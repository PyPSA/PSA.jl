using JuMP

include("auxilliaries.jl")

function lopf(network, solver, formulation, objective, investment)
    # This function is organized as the following:
    #
    # 0.        Initialize model
    # 1. - 5.   add generators,  lines, links, storage_units,
    #           stores to the model:
    #               .1 separate different types from each other
    #               .2 define number of different types
    #               .3 add variables to the model
    #               .4 set contraints for extendables
    #               .5 set charging constraints (storage_units and stores)
    # 6.        set Kirchhoff's Current Law (nodal balance)
    # 7.        set Kirchhoff Voltage Law
    # 8.        set global constraints
    # 9.        objective function and solving
    # 10.       extract results

    # Conventions:
    #   - all variable names start with capital letters
    #   - Additional capital words are N (number of buses), T (number of snapshots)
    #       and N_* defining the nuber of considered variables added to the model.

    #solver is e.g.

    #using Gurobi
    #solver = GurobiSolver(Method=1,Threads=2)

    #using Clp
    #solver = ClpSolver()

    println("Creating model.")
    m = Model(solver=solver)

    calculate_dependent_values!(network)
    buses = network.buses
    lines = network.lines
    reverse_busidx = rev_idx(buses)
    busidx = idx(buses)
    reverse_lineidx = rev_idx(lines)
    lineidx = idx(lines)
    N = nrow(network.buses)
    L = nrow(network.lines)
    T = nrow(network.snapshots) #normally snapshots
    nrow(network.loads_t["p"])!=T ? network.loads_t["p"]=network.loads_t["p_set"] : nothing

# --------------------------------------------------------------------------------------------------------

# 1. add all generators to the model
    println("Adding generators to the model.")

    # 1.1 set different generator types
    println("-- 1.1")
    generators = network.generators
    fix_gens_b = ((.!generators[:p_nom_extendable]) .& (.!generators[:commitable]))
    ext_gens_b = convert(BitArray, generators[:p_nom_extendable])
    com_gens_b = convert(BitArray, generators[:commitable])

    # 1.2 fix bounds for iterating
    println("-- 1.2")
    N_fix = sum(fix_gens_b)
    N_ext = sum(ext_gens_b)
    N_com = sum(com_gens_b)

    if N_com > 0
        println("WARNING, no unit commitment yet")
    end

    p_max_pu = get_switchable_as_dense(network, "generators", "p_max_pu")
    p_min_pu = get_switchable_as_dense(network, "generators", "p_min_pu")

    # 1.3a add non-extendable generator variables to the model
    println("-- 1.3a")
    p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=fix_gens_b)
    p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=fix_gens_b)
    p_nom = network.generators[fix_gens_b,:p_nom]

    @variable m p_nom[gr]*p_min_pu(t,gr) <= G_fix[gr=1:N_fix,t=1:T] <= p_nom[gr]*p_max_pu(t,gr)


    # 1.3b add extendable generator variables to the model
    println("-- 1.3b")
    p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=ext_gens_b)
    p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=ext_gens_b)
    p_nom_min = network.generators[ext_gens_b,:p_nom_min]
    p_nom_max = network.generators[ext_gens_b,:p_nom_max]

    @variable m G_ext[gr=1:N_ext,t = 1:T]

    @variable m p_nom_min[gr] <=  gen_p_nom[gr=1:N_ext] <= p_nom_max[gr]

    @constraints m begin
        [gr=1:N_ext,t=1:T], G_ext[gr,t] >= gen_p_nom[gr]*p_min_pu(t,gr)
        [gr=1:N_ext,t=1:T], G_ext[gr,t] <= gen_p_nom[gr]*p_max_pu(t,gr)
    end


    G = [G_fix; G_ext] # G is the concatenated variable array
    generators = [generators[fix_gens_b,:]; generators[ext_gens_b,:] ] # sort generators the same
    # new booleans
    fix_gens_b = ((.!generators[:p_nom_extendable]) .& (.!generators[:commitable]))
    ext_gens_b = convert(BitArray, generators[:p_nom_extendable])
    com_gens_b = convert(BitArray, generators[:commitable])



# --------------------------------------------------------------------------------------------------------

# 2. add all lines to the model
    println("Adding lines to the model.")

    # 2.1 set different lines types
    println("-- 2.1")
    lines = network.lines
    fix_lines_b = (.!lines[:s_nom_extendable])
    ext_lines_b = .!fix_lines_b

    # 2.2 iterator bounds
    println("-- 2.2")
    N_fix = sum(fix_lines_b)
    N_ext = sum(ext_lines_b)

    # 2.3 add variables
    println("-- 2.3")
    @variables m begin
        -lines[fix_lines_b,:s_nom][l]  <=  LN_fix[l=1:N_fix,t=1:T] <= lines[fix_lines_b,:s_nom][l]
        LN_ext[l=1:N_ext,t=1:T]
        lines[ext_lines_b,:s_nom_min][l] <=  LN_s_nom[l=1:N_ext] <= lines[ext_lines_b,:s_nom_max][l]
    end

    # 2.4 add line constraint for extendable lines
    println("-- 2.4")
    @constraints(m, begin
            [l=1:N_ext,t=1:T], LN_ext[l,t] <=  LN_s_nom[l]
            [l=1:N_ext,t=1:T], LN_ext[l,t] >= -LN_s_nom[l]
    end)

    # 2.5 add integer variables
    if investment == "int"
        @variable(m, LN_opt_inv[l=1:N_ext], Int)
        @constraint(m, integer[l=1:N_ext], LN_s_nom[l] == LN_opt_inv[l] * lines[ext_lines_b,:s_nom_step][l] + lines[ext_lines_b,:s_nom][l])

    elseif investment == "bin_abs"
        bigM = min.(lines[ext_lines_b,:s_nom_max],1e6) # TODO choose bigM default
        @variable(m, LN_opt[l=1:N_ext], Bin)
        @variable(m, LN_inv[l=1:N_ext])
        @constraints m begin
            [l=1:N_ext], - bigM[l] * (1-LN_opt[l]) + lines[ext_lines_b,:abs_ext_min][l] <= LN_inv[l]
            [l=1:N_ext], - bigM[l] * LN_opt[l] <= LN_inv[l]
            [l=1:N_ext],  bigM[l] * LN_opt[l] >= LN_inv[l]
            [l=1:N_ext], LN_s_nom[l] == LN_inv[l] + lines[ext_lines_b,:s_nom][l]
        end

    elseif investment == "convexhull"
        bigM = min.(lines[ext_lines_b,:s_nom_max],1e6) # TODO choose bigM default
        @variable(m, LN_opt[l=1:N_ext], Bin)
        @variable(m, LN_inv[l=1:N_ext])
        @variable(m, LN_inv_1[l=1:N_ext])
        @variable(m, LN_inv_2[l=1:N_ext])
        #@constraint(m, [l=1:N_ext], LN_inv[l] = LN_inv_1[l] + LN_inv_2[l])
        @constraints m begin
            [l=1:N_ext], 0 <= LN_inv_1[l]
            [l=1:N_ext], LN_inv_1[l] <= bigM[l] * (1-LN_opt[l])
            [l=1:N_ext], 0 <= LN_inv_2[l]
            [l=1:N_ext], LN_inv_2[l] <= bigM[l] * LN_opt[l]
        end
        @constraints m begin
            [l=1:N_ext], (1-LN_opt[l]) * lines[ext_lines_b,:abs_ext_min][l] <= LN_inv_1[l]
            [l=1:N_ext], 0 == LN_inv_2[l]
            [l=1:N_ext], LN_s_nom[l] == LN_inv_1[l] + LN_inv_2[l] + lines[ext_lines_b,:s_nom][l]
        end

    # TODO requires non infinite s_nom_max values to work properly (and make sense):
    # current workaround: choose very high s_nom_max as default and set rel_ext_min to 0.0
    # this way minimum expansion constraint is ignored
    elseif investment == "bin_rel"

        # TODO catch infinite s_nom_max lines, assigning new values in DataFrame not working yet
        s_nom_max = lines[ext_lines_b,:s_nom_max]
        rel_ext_min = lines[ext_lines_b,:rel_ext_min]
        for l=1:N_ext
            if s_nom_max[l] == Inf
                s_nom_max[l] = 1e6 # TODO not working ugly workaround
                rel_ext_min[l] = 0.0
            end
        end

        bigM = min.(lines[ext_lines_b,:s_nom_max],1e6) # TODO choose bigM default
        @variable(m, LN_opt[l=1:N_ext], Bin)
        @variable(m, LN_inv[l=1:N_ext])
        @constraints m begin
            [l=1:N_ext], - bigM[l] * (1-LN_opt[l]) * rel_ext_min[l] <= LN_inv[l]
            [l=1:N_ext], - bigM[l] * LN_opt[l] <= LN_inv[l]
            [l=1:N_ext],  bigM[l] * LN_opt[l] >= LN_inv[l]
            [l=1:N_ext], LN_s_nom[l] == LN_inv[l]* ( s_nom_max[l]-lines[ext_lines_b,:s_nom][l] ) + lines[ext_lines_b,:s_nom][l]
        end

    end

    LN = [LN_fix; LN_ext]
    lines = [lines[fix_lines_b,:]; lines[ext_lines_b,:]]
    fix_lines_b = (.!lines[:s_nom_extendable])
    ext_lines_b = .!fix_lines_b

# --------------------------------------------------------------------------------------------------------

# 3. add all links to the model
    println("Adding links to the model.")

    # 3.1 set different link types
    println("-- 3.1")
    links = network.links
    fix_links_b = .!links[:p_nom_extendable]
    ext_links_b = .!fix_links_b

    # 3.2 iterator bounds
    println("-- 3.2")
    N_fix = sum(fix_links_b)
    N_ext = sum(ext_links_b)

    #  3.3 set link variables
    println("-- 3.3")
    @variables m begin
       ((links[fix_links_b, :p_nom].*links[fix_links_b, :p_min_pu])[l]  <=  LK_fix[l=1:N_fix,t=1:T]
                <= (links[fix_links_b, :p_nom].*links[fix_links_b, :p_max_pu])[l])
        LK_ext[l=1:N_ext,t=1:T]
        links[ext_links_b, :p_nom_min][l] <=  LK_p_nom[l=1:N_ext] <= links[ext_links_b, :p_nom_max][l]
    end

    # 3.4 set constraints for extendable links
    println("-- 3.4")
    @constraints(m, begin
    [l=1:N_ext,t=1:T], LK_ext[l,t] >= LK_p_nom[l].*links[ext_links_b, :p_min_pu][l]
    [l=1:N_ext,t=1:T], LK_ext[l,t] <= LK_p_nom[l].*links[ext_links_b, :p_max_pu][l]
    end)

    LK = [LK_fix; LK_ext]
    links = [links[fix_links_b,:]; links[ext_links_b,:]]
    fix_links_b = .!links[:p_nom_extendable]
    ext_links_b = .!fix_links_b

# --------------------------------------------------------------------------------------------------------
# 4. define storage_units
    println("Adding storage units to the model.")

    # 4.1 set different storage_units types
    println("-- 4.1")
    storage_units = network.storage_units
    fix_sus_b = .!storage_units[:p_nom_extendable]
    ext_sus_b = .!fix_sus_b

    inflow = get_switchable_as_dense(network, "storage_units", "inflow")

    # 4.2 iterator bounds
    println("-- 4.2")
    N_fix = sum(fix_sus_b)
    N_ext = sum(ext_sus_b)
    N_sus = nrow(storage_units)

    #  4.3 set variables
    println("-- 4.3")
    @variables m begin
       (0 <=  SU_dispatch_fix[s=1:N_fix,t=1:T] <=
                (storage_units[fix_sus_b, :p_nom].*storage_units[fix_sus_b, :p_max_pu])[s])

        SU_dispatch_ext[s=1:N_ext,t=1:T] >= 0

        (0 <=  SU_store_fix[s=1:N_fix,t=1:T] <=
                 - (storage_units[fix_sus_b, :p_nom].*storage_units[fix_sus_b, :p_min_pu])[s])

        SU_store_ext[s=1:N_ext,t=1:T] >= 0

        SU_p_nom[s=1:N_ext] >= 0

        0 <= SU_soc_fix[s=1:N_fix,t=1:T] <= (storage_units[fix_sus_b,:max_hours] #  TODO [s]?
                                            .*storage_units[fix_sus_b,:p_nom])[s]

        SU_soc_ext[s=1:N_ext,t=1:T] >= 0

        0 <=  SU_spill_fix[s=1:N_fix,t=1:T] <= inflow[:,fix_sus_b][t,s]

        0 <=  SU_spill_ext[s=1:N_ext,t=1:T] <= inflow[:,ext_sus_b][t,s]
    end

    # 4.4 set constraints for extendable storage_units
    println("-- 4.4")

    @constraints(m, begin
            [s=1:N_ext,t=1:T], SU_dispatch_ext[s,t] <= SU_p_nom[s].*storage_units[ext_sus_b, :p_max_pu][s]
            [s=1:N_ext,t=1:T], SU_store_ext[s,t] <= - SU_p_nom[s].*storage_units[ext_sus_b, :p_min_pu][s]
            [s=1:N_ext,t=1:T], SU_soc_ext[s,t] <= SU_p_nom[s].*storage_units[ext_sus_b, :max_hours][s]
    end)

    # 4.5 set charging constraint
    println("-- 4.5")
    SU_dispatch = [SU_dispatch_fix; SU_dispatch_ext]
    SU_store = [SU_store_fix; SU_store_ext]
    SU_soc = [SU_soc_fix; SU_soc_ext]
    SU_spill = [SU_spill_fix; SU_spill_ext]

    storage_units = [storage_units[fix_sus_b,:]; storage_units[ext_sus_b,:]]
    inflow = [inflow[:,fix_sus_b] inflow[:,ext_sus_b]]

    ext_sus_b = BitArray(storage_units[:p_nom_extendable])

    is_cyclic_i = collect(1:N_sus)[BitArray(storage_units[:cyclic_state_of_charge])]
    not_cyclic_i = collect(1:N_sus)[.!storage_units[:cyclic_state_of_charge]]

    @constraints(m, begin
            [s=is_cyclic_i], SU_soc[s,1] == (SU_soc[s,T]
                                        + storage_units[s,:efficiency_store] * SU_store[s,1]
                                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,1]
                                        + inflow[1,s] - SU_spill[s,1] )
            [s=not_cyclic_i], SU_soc[s,1] == (storage_units[s,:state_of_charge_initial]
                                        + storage_units[s,:efficiency_store] * SU_store[s,1]
                                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,1]
                                        + inflow[1,s] - SU_spill[s,1])

            [s=1:N_sus,t=2:T], SU_soc[s,t] == (SU_soc[s,t-1]
                                            + storage_units[s,:efficiency_store] * SU_store[s,t]
                                            - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,t]
                                            + inflow[t,s] - SU_spill[s,t] )

        end)

# --------------------------------------------------------------------------------------------------------

# 5. define stores
    println("Adding stores to the model.")

# 5.1 set different stores types
    println("-- 5.1")
    stores = network.stores

    fix_stores_b = .!stores[:e_nom_extendable]
    ext_stores_b = .!fix_stores_b

    inflow = get_switchable_as_dense(network, "stores", "inflow")

    # 5.2 iterator bounds
    println("-- 5.2")
    N_fix = sum(fix_stores_b)
    N_ext = sum(ext_stores_b)
    N_st = N_fix + N_ext


    #  5.3 set variables
    println("-- 5.3")
    @variables m begin
       (0 <=  ST_dispatch_fix[s=1:N_fix,t=1:T] <=
                (stores[fix_stores_b, :e_nom].*stores[fix_stores_b, :e_max_pu])[s])
        ST_dispatch_ext[s=1:N_ext,t=1:T] >= 0
        (0 <=  ST_store_fix[s=1:N_fix,t=1:T] <=
                 - (stores[fix_stores_b, :e_nom].*stores[fix_stores_b, :e_min_pu])[s])
        ST_store_ext[s=1:N_ext,t=1:T] >= 0

        ST_e_nom[s=1:N_ext] >= 0

        0 <= ST_soc_fix[s=1:N_fix,t=1:T] <= (stores[fix_stores_b,:max_hours]
                                            .*stores[fix_stores_b,:e_nom])[s]
        ST_soc_ext[s=1:N_ext,t=1:T] >= 0

        0 <=  ST_spill_fix[l=1:N_fix,t=1:T] <= inflow[:,fix_stores_b][l=1:N_fix,t=1:T]
        0 <=  ST_spill_ext[l=1:N_fix,t=1:T] <= inflow[:,ext_stores_b][l=1:N_ext,t=1:T]
    end

    # 5.4 set constraints for extendable stores
    println("-- 5.4")

    @constraints(m, begin
            [s=1:N_ext,t=1:T], ST_dispatch_ext[s,t] <= ST_e_nom[s].*stores[ext_stores_b, :e_max_pu][s]
            [s=1:N_ext,t=1:T], ST_store_ext[s,t] <= - ST_e_nom[s].*stores[ext_stores_b, :e_min_pu][s]
            [s=1:N_ext,t=1:T], ST_soc_ext[s,t] <= ST_e_nom[s].*stores[ext_stores_b, :max_hours][s]
    end)

    # 5.5 set charging constraint
    println("-- 5.5")
    ST_dispatch = [ST_dispatch_fix; ST_dispatch_ext]
    ST_store = [ST_store_fix; ST_store_ext]
    ST_soc = [ST_soc_fix; ST_soc_ext]
    ST_spill = [ST_spill_fix, ST_spill_ext]
    stores = [stores[fix_stores_b,:]; stores[ext_stores_b,:]]

    inflow = [inflow[:,fix_stores_b] inflow[:,ext_stores_b]]

    ext_stores_b = BitArray(stores[:e_nom_extendable])

    is_cyclic_i = collect(1:N_st)[BitArray(stores[:cyclic_state_of_charge])]
    not_cyclic_i = collect(1:N_st)[.!stores[:cyclic_state_of_charge]]

    @constraints(m, begin
            [s=is_cyclic_i,t=1], ST_soc[s,t] == (ST_soc[s,T]
                                        + stores[s,:efficiency_store] * ST_store[s,t]
                                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t] # TODO this is different to SU
                                        + inflow[t,s] - ST_spill[s,t])

            [s=not_cyclic_i,t=1], ST_soc[s,t] == (stores[s,:state_of_charge_initial]
                                        + stores[s,:efficiency_store] * ST_store[s,t]
                                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]  # TODO this is different to SU
                                        + inflow[t,s] - ST_spill[s,t])

            [s=is_cyclic_i,t=2:T], ST_soc[s,t] == (ST_soc[s,t-1]
                                            + stores[s,:efficiency_store] * ST_store[s,t]
                                            - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]  # TODO this is different to SU
                                            + inflow[t,s] - ST_spill[s,t])

        end)

# --------------------------------------------------------------------------------------------------------
# 6. + 7. power flow formulations
    # TODO need validation / testing (e.g. why is there a difference for 37 node example; why does kirchhoff formulation deviate?)

    println("Adding power flow formulation $formulation to the model.")

    # a. angles formulation
    if formulation == "angles"

        # voltage angles
        @variable(m, theta[1:N,1:T])

        K = incidence_matrix(network)
        Lambda = weighted_laplace_matrix(network)
        B = reactance_matrix(network)

        # load data in correct order
        loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

        @constraint(m, voltages[n=1:N, t=1:T], (

                  sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                      .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                == Lambda[n,:]' * theta[:,t] ))

        @constraint(m, flows[l=1:L, t=1:T], LN[l, t] == (B * K' * theta[:,t])[l] )

        @constraint(m, slack[t=1:T], theta[1,t] == 0 )

    # b. cycles formulation
    elseif formulation == "cycles"

        #load data in correct order
        loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

        @constraint(m, balance[n=1:N, t=1:T], (
              sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
            + sum(LN[ findin(lines[:bus1], [reverse_busidx[n]]) ,t])
            + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                  .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
            + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

            - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
            - sum(LN[ findin(lines[:bus0], [reverse_busidx[n]]) ,t])
            - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
            - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

              == 0 ))

        # since cyclebasis is not yet supported in LightGraphs, use the pyimported
        # netwrokx in order to define all cycles. The cycle_basis returns a list of
        # cycles, indicating the connected buses. For each cycle the connecting branches
        # and their directions (signs) have to be determined.

        # Might be nessecary to loop over all subgraphs as
        # for (sn, sub) in enumerate(weakly_connected_components(g))
        #     g_sub = induced_subgraph(g, sub)[1]

        append_idx_col!(lines)
        (branches, var, attribute) = (lines, LN, :x)
        cycles = get_cycles(network)
        if ndims(cycles)<2
            cycles = [cycle for cycle in cycles if length(cycle)>2]
        else
            cycles = [cycles[i,:] for i in 1:size(cycles)[1]]
        end
        if length(cycles)>0
            cycles_branch = Array{Int64,1}[]
            directions = Array{Float64,1}[]
            for cyc=1:length(cycles)
                push!(cycles_branch,Int64[])
                push!(directions,Float64[])
                for bus=1:length(cycles[cyc])
                    bus0 = cycles[cyc][bus]
                    bus1 = cycles[cyc][(bus)%length(cycles[cyc])+1]
                    try
                        push!(cycles_branch[cyc],branches[((branches[:bus0].==reverse_busidx[bus0])
                                    .&(branches[:bus1].==reverse_busidx[bus1])),:idx][1] )
                        push!(directions[cyc], 1.)
                    catch y
                        if isa(y, BoundsError)
                            push!(cycles_branch[cyc], branches[((branches[:bus0].==reverse_busidx[bus1])
                                            .&(branches[:bus1].==reverse_busidx[bus0])),:idx][1] )
                            push!(directions[cyc], -1.)
                        else
                            return y
                        end
                    end
                end
            end
            if attribute==:x
                @constraint(m, line_cycle_constraint[c=1:length(cycles_branch), t=1:T] ,
                        dot(directions[c] .* lines[cycles_branch[c], :x_pu],
                            LN[cycles_branch[c],t]) == 0)
            # elseif attribute==:r
            #     @constraint(m, link_cycle_constraint[c=1:length(cycles_branch), t=1:T] ,
            #             dot(directions[c] .* links[cycles_branch[c], :r]/380. , LK[cycles_branch[c],t]) == 0)
            end
        end

    # c. kirchhoff formulation
    elseif formulation == "kirchhoff"

        K = incidence_matrix(network)
        C = get_cycles(network)

        # adapting C to handle 1 cycle networks
        ndims(C) == 2 ? C = [C] : nothing

        n_cycles = length(C)
        cycles = zeros(Int64, L, n_cycles)

        for c=1:n_cycles
            for l in C[c]
                cycles[l,c] = 1
            end
        end

        #load data in correct order
        loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

        @constraint(m, kcl[n=1:N,t=1:T],
            (    sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
               + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                 .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
               + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

               - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
               - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
               - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
             ) == K[n,:]' * LN[:,t] )

        @constraint(m, kvl[c=1:n_cycles,t=1:T],
            sum(cycles[l,c] * network.lines[:x][l] * LN[l,t] for l=1:L) == 0 )

    # d. ptdf formulation
    elseif formulation == "ptdf"

        ptdf = ptdf_matrix(network)

        #load data in correct order
        loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

        @constraint(m, flows[l=1:L,t=1:T],
            sum( ptdf[l,n]
               * (   sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                   + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                     .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                   + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                   - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                   - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                   - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
                 )
               for n in 1:N
               ) == LN[l,t] )

        @constraint(m, balance[t=1:T],
        sum( (   sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
               + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                 .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
               + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

               - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
               - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
               - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
             ) for n in 1:N
           ) == 0 )

    else
        println("The formulation $formulation is not implemented.")
    end

# --------------------------------------------------------------------------------------------------------

# 8. set global_constraints
# only for co2_emissions till now
    println("Adding global constraints to the model.")

    if nrow(network.global_constraints)>0 && in("primary_energy", network.global_constraints[:type])
        co2_limit = network.global_constraints[network.global_constraints[:name].=="co2_limit", :constant]
        nonnull_carriers = network.carriers[network.carriers[:co2_emissions].!=0, :]
        # emmssions = Dict(zip(nonnull_carriers[:name], nonnull_carriers[:co2_emissions])) // unused
        carrier_index(carrier) = findin(generators[:carrier], [carrier])
        @constraint(m, sum(sum(dot(1./generators[carrier_index(carrier) , :efficiency],
                    G[carrier_index(carrier),t]) for t=1:T)
                    * select_names(network.carriers, [carrier])[:co2_emissions]
                    for carrier in network.carriers[:name]) .<=  co2_limit)
    end

# --------------------------------------------------------------------------------------------------------

# 9. set objective function
    println("Adding objective to the model.")

    # TODO catch if some components are not existing in network!
    # TODO this objective consider the total cost (capex of existing + capex of additional) of extensible components!
    # Should be only or extensions.
    # Accurate as long as extensible components have 0 capacity at beginning!
    if objective == "total"
        @objective(m, Min,
                              sum(dot(generators[:marginal_cost], G[:,t]) for t=1:T)
                            + dot(generators[ext_gens_b,:capital_cost], gen_p_nom[:] )
                            + dot(generators[fix_gens_b,:capital_cost], generators[fix_gens_b,:p_nom])

                            + dot(lines[ext_lines_b,:capital_cost], LN_s_nom[:])
                            + dot(lines[fix_lines_b,:capital_cost], lines[fix_lines_b,:s_nom])

                            + dot(links[ext_links_b,:capital_cost], LK_p_nom[:])
                            + dot(links[fix_links_b,:capital_cost], links[fix_links_b,:p_nom])

                            + sum(dot(storage_units[:marginal_cost], SU_dispatch[:,t]) for t=1:T)
                            + dot(storage_units[ext_sus_b, :capital_cost], SU_p_nom[:])
                            + dot(storage_units[fix_sus_b,:capital_cost], storage_units[fix_sus_b,:p_nom])

                            + sum(dot(stores[:marginal_cost], ST_dispatch[:,t]) for t=1:T)
                            + dot(stores[ext_stores_b, :capital_cost], ST_e_nom[:])
                            + dot(stores[fix_stores_b,:capital_cost], stores[fix_stores_b,:e_nom])
                    )
    elseif objective == "extensions"
        # not working yet!
        # TODO modify objective function such that only cost of extensions and opex is included
        # reducing capacity saves money; in reality would increase cost
        @objective(m, Min,
                              sum(dot(generators[:marginal_cost], G[:,t]) for t=1:T)
                            + dot(generators[ext_gens_b,:capital_cost], gen_p_nom[:] - generators[ext_gens_b,:p_nom])

                            + dot(lines[ext_lines_b,:capital_cost], LN_s_nom[:] - lines[ext_lines_b,:s_nom])

                            + dot(links[ext_links_b,:capital_cost], LK_p_nom[:] - links[ext_links_b,:p_nom])

                            + sum(dot(storage_units[:marginal_cost], SU_dispatch[:,t]) for t=1:T)
                            + dot(storage_units[ext_sus_b, :capital_cost], SU_p_nom[:] - storage_units[ext_sus_b,:p_nom])

                            + sum(dot(stores[:marginal_cost], ST_dispatch[:,t]) for t=1:T)
                            + dot(stores[ext_stores_b, :capital_cost], ST_e_nom[:] - stores[ext_stores_b,:e_nom])
                    )
    else # TODO original objective; delete when other objectives tested.
        @objective(m, Min,
                              sum(dot(generators[:marginal_cost], G[:,t]) for t=1:T)
                            + dot(generators[ext_gens_b,:capital_cost], gen_p_nom[:]  )

                            + dot(lines[ext_lines_b,:capital_cost], LN_s_nom[:])

                            + dot(links[ext_links_b,:capital_cost], LK_p_nom[:])

                            + sum(dot(storage_units[:marginal_cost], SU_dispatch[:,t]) for t=1:T)
                            + dot(storage_units[ext_sus_b, :capital_cost], SU_p_nom[:])

                            + sum(dot(stores[:marginal_cost], ST_dispatch[:,t]) for t=1:T)
                            + dot(stores[ext_stores_b, :capital_cost], ST_e_nom[:])

                    )
    end


    status = solve(m)

# --------------------------------------------------------------------------------------------------------

# 10. extract optimisation results
    if status==:Optimal

        orig_gen_order = network.generators[:name]
        generators[:p_nom_opt] = deepcopy(generators[:p_nom])
        generators[ext_gens_b,:p_nom_opt] = getvalue(gen_p_nom)
        network.generators = generators
        network.generators_t["p"] = names!(DataFrame(transpose(getvalue(G))), Symbol.(generators[:name]))
        network.generators = select_names(network.generators, orig_gen_order)

        orig_line_order = network.lines[:name]
        network.lines = lines
        lines[:s_nom_opt] = deepcopy(lines[:s_nom])
        network.lines[ext_lines_b,:s_nom_opt] = getvalue(LN_s_nom)
        network.lines_t["p0"] = names!(DataFrame(transpose(getvalue(LN))), Symbol.(lines[:name]))
        network.lines = select_names(network.lines, orig_line_order)

        # network.buses_t["p"] =  DataFrame(ncols=nrow(network.buses))

        if nrow(links)>0
            orig_link_order = network.links[:name]
            network.links = links
            links[:p_nom_opt] = deepcopy(links[:p_nom])
            network.links[ext_links_b,:p_nom_opt] = getvalue(LK_p_nom)
            network.links_t["p0"] = names!(DataFrame(transpose(getvalue(LK))), Symbol.(links[:name]))
            network.links = select_names(network.links, orig_link_order)

        end
        if nrow(storage_units)>0
            orig_sus_order = network.storage_units[:name]
            network.storage_units = storage_units

            storage_units[:p_nom_opt] = deepcopy(storage_units[:p_nom])
            network.storage_units[ext_sus_b,:p_nom_opt] = getvalue(SU_p_nom)
            # network.storage_units_t["spill"] = spillage
            # network.storage_units_t["spill"][:,spill_sus_b] = names!(DataFrame(transpose(getvalue(SU_spill))),
            #                     names(spillage)[spill_sus_b])
            network.storage_units_t["spill"] = names!(DataFrame(transpose(getvalue(SU_spill))),
                                Symbol.(storage_units[:name]))
            network.storage_units_t["p"] = names!(DataFrame(transpose(getvalue(SU_dispatch .- SU_store))),
                                Symbol.(storage_units[:name]))
            network.storage_units_t["state_of_charge"] = names!(DataFrame(transpose(getvalue(SU_soc))),
                                Symbol.(storage_units[:name]))
            network.storage_units = select_names(network.storage_units, orig_sus_order)
        end
        align_component_order!(network)
        println("Reduce cost to $(m.objVal)")
    end
    return m
end
