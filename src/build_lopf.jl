using JuMP
using MathProgBase

include("utils.jl")

function build_lopf(network, solver; formulation::String="angles", objective::String="total", investment_type::String="continuous", blockmodel::Bool=false)
    
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
    # 6.        give power flow formulation
    # 7.        set global constraints
    # 8.        give objective function

    # Conventions:
    #   - all variable names start with capital letters
    #   - Additional capital words - which are not variables - are N (number of buses),
    #     T (number of snapshots), L (number of lines) defining the number of considered variables added to the model.

    println("Creating model.")

    if blockmodel
        m = BlockModel(solver=solver)
    else
        m = Model(solver=solver)
    end

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
    #println("-- 1.1 set different generator types")
    generators = network.generators
    fix_gens_b = ((.!generators[:p_nom_extendable]) .& (.!generators[:commitable]))
    ext_gens_b = convert(BitArray, generators[:p_nom_extendable])
    com_gens_b = convert(BitArray, generators[:commitable])

    # 1.2 fix bounds for iterating
    #println("-- 1.2 fix bounds for iterating")
    N_fix = sum(fix_gens_b)
    N_ext = sum(ext_gens_b)
    N_com = sum(com_gens_b)

    if N_com > 0
        println("WARNING, no unit commitment yet")
    end

    p_max_pu = get_switchable_as_dense(network, "generators", "p_max_pu")
    p_min_pu = get_switchable_as_dense(network, "generators", "p_min_pu")

    # 1.3a add non-extendable generator variables to the model
    #println("-- 1.3a add non-extendable generator variables to the model")
    p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=fix_gens_b)
    p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=fix_gens_b)
    p_nom = network.generators[fix_gens_b,:p_nom]

    @variable m p_nom[gr]*p_min_pu(t,gr) <= G_fix[gr=1:N_fix,t=1:T] <= p_nom[gr]*p_max_pu(t,gr)

    # 1.3b add extendable generator variables to the model
    #println("-- 1.3b add extendable generator variables to the model")
    p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=ext_gens_b)
    p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=ext_gens_b)
    p_nom_min = network.generators[ext_gens_b,:p_nom_min]
    p_nom_max = network.generators[ext_gens_b,:p_nom_max]

    @variable m G_ext[gr=1:N_ext,t = 1:T]

    @variable m p_nom_min[gr] <=  G_p_nom[gr=1:N_ext] <= p_nom_max[gr]

    @constraints m begin
        lower_gen_limit[gr=1:N_ext,t=1:T], G_ext[gr,t] >= G_p_nom[gr]*p_min_pu(t,gr)
        upper_gen_limit[gr=1:N_ext,t=1:T], G_ext[gr,t] <= G_p_nom[gr]*p_max_pu(t,gr)
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
    #println("-- 2.1 set different lines types")
    lines = network.lines
    fix_lines_b = (.!lines[:s_nom_extendable])
    ext_lines_b = .!fix_lines_b

    # 2.2 iterator bounds
    #println("-- 2.2 iterator bounds")
    N_fix = sum(fix_lines_b)
    N_ext = sum(ext_lines_b)

    # 2.3 add variables
    #println("-- 2.3 add variables")
    @variables m begin
        -lines[fix_lines_b,:s_nom][l]  <=  LN_fix[l=1:N_fix,t=1:T] <= lines[fix_lines_b,:s_nom][l]
        LN_ext[l=1:N_ext,t=1:T]
        lines[ext_lines_b,:s_nom_min][l] <=  LN_s_nom[l=1:N_ext] <= lines[ext_lines_b,:s_nom_max][l]
    end

    # 2.4 add line constraint for extendable lines
    #println("-- 2.4 add line constraint for extendable lines")
    @constraints(m, begin
            upper_line_limit[l=1:N_ext,t=1:T], LN_ext[l,t] <=  LN_s_nom[l]
            lower_line_limit[l=1:N_ext,t=1:T], LN_ext[l,t] >= -LN_s_nom[l]
    end)

    # 2.5 add integer variables if applicable
    #println("-- 2.5 add integer variables if applicable")

    bigM_default = 1e4 # TODO: choose bigM default

    if investment_type == "continuous"
        @variable(m, LN_inv[l=1:N_ext]) 
        #@constraint(m, continuous[l=1:N_ext], LN_s_nom[l] == LN_inv[l] * lines[ext_lines_b,:s_nom_step][l] + lines[ext_lines_b,:s_nom][l]) # if s_nom_step is specified and option "angles_bilinear_stepsspecs"
        @constraint(m, continuous[l=1:N_ext], LN_s_nom[l] == (1+LN_inv[l]) * lines[ext_lines_b,:s_nom][l])

    elseif investment_type == "integer"
        @variable(m, LN_inv[l=1:N_ext], Int) 
        #@constraint(m, integer[l=1:N_ext], LN_s_nom[l] == LN_inv[l] * lines[ext_lines_b,:s_nom_step][l] + lines[ext_lines_b,:s_nom][l]) # if s_nom_step is specified and option "angles_bilinear_stepsspecs
        @constraint(m, integer[l=1:N_ext], LN_s_nom[l] == (1+LN_inv[l]) * lines[ext_lines_b,:s_nom][l])
    
    # if investment in line is chosen, it must be above a minimum threshold s_nom_ext_min
    # apart from that investment is continuous
    elseif investment_type == "binary"
        bigM = min.(lines[ext_lines_b,:s_nom_max],bigM_default)
        @variable(m, LN_opt[l=1:N_ext], Bin)
        @variable(m, LN_inv[l=1:N_ext])
        @constraints(m, begin
            binary1[l=1:N_ext], - bigM[l] * (1-LN_opt[l]) + lines[ext_lines_b,:s_nom_ext_min][l] <= LN_inv[l]
            binary2[l=1:N_ext], 0 <= LN_inv[l] # no reduction of capacity allowed
            binary3[l=1:N_ext],  bigM[l] * LN_opt[l] >= LN_inv[l]
            binary4[l=1:N_ext], LN_s_nom[l] == (1+LN_inv[l]) * lines[ext_lines_b,:s_nom][l]
        end)

    elseif investment_type == "integer_bigm"
        
        candidates = Array{Int64,1}[]
        for l=1:N_ext
            if lines[:s_nom_max][l] != Inf
                max_extension = floor(lines[:s_nom_max][l] / lines[:s_nom][l])
                push!(candidates,[i for i=1:max_extension])
            else
                # fallback option
                push!(candidates,[i for i=1:10])
            end
        end
        println(candidates)

        @variable(m, LN_opt[l=1:N_ext,c in candidates[l]], Bin)

        @constraint(m, logical[l=1:N_ext], sum(LN_opt[l,c] for c in candidates[l]) == 1)

        @constraint(m, integer_binary_reformulation[l=1:N_ext], LN_s_nom[l] == (sum(c*LN_opt[l,c] for c in candidates[l])) * lines[ext_lines_b,:s_nom][l])

    end

    LN = [LN_fix; LN_ext]
    lines = [lines[fix_lines_b,:]; lines[ext_lines_b,:]]
    fix_lines_b = (.!lines[:s_nom_extendable])
    ext_lines_b = .!fix_lines_b

# --------------------------------------------------------------------------------------------------------

# 3. add all links to the model
    println("Adding links to the model.")

    # 3.1 set different link types
    #println("-- 3.1 set different link types")
    links = network.links
    fix_links_b = .!links[:p_nom_extendable]
    ext_links_b = .!fix_links_b

    # 3.2 iterator bounds
    #println("-- 3.2 iterator bounds")
    N_fix = sum(fix_links_b)
    N_ext = sum(ext_links_b)

    #  3.3 set link variables
    #println("-- 3.3 set link variables")
    @variables m begin
       ((links[fix_links_b, :p_nom].*links[fix_links_b, :p_min_pu])[l]  <=  LK_fix[l=1:N_fix,t=1:T]
                <= (links[fix_links_b, :p_nom].*links[fix_links_b, :p_max_pu])[l])
        LK_ext[l=1:N_ext,t=1:T]
        links[ext_links_b, :p_nom_min][l] <=  LK_p_nom[l=1:N_ext] <= links[ext_links_b, :p_nom_max][l]
    end

    # 3.4 set constraints for extendable links
    #println("-- 3.4 set constraints for extendable links")
    @constraints(m, begin
        lower_link_limit[l=1:N_ext,t=1:T], LK_ext[l,t] >= LK_p_nom[l].*links[ext_links_b, :p_min_pu][l]
        upper_link_limit[l=1:N_ext,t=1:T], LK_ext[l,t] <= LK_p_nom[l].*links[ext_links_b, :p_max_pu][l]
    end)

    LK = [LK_fix; LK_ext]
    links = [links[fix_links_b,:]; links[ext_links_b,:]]
    fix_links_b = .!links[:p_nom_extendable]
    ext_links_b = .!fix_links_b

# --------------------------------------------------------------------------------------------------------

# 4. define storage_units
    println("Adding storage units to the model.")

    # 4.1 set different storage_units types
    #println("-- 4.1 set different storage_units types")
    storage_units = network.storage_units
    fix_sus_b = .!storage_units[:p_nom_extendable]
    ext_sus_b = .!fix_sus_b

    inflow = get_switchable_as_dense(network, "storage_units", "inflow")

    # 4.2 iterator bounds
    #println("-- 4.2 iterator bounds")
    N_fix = sum(fix_sus_b)
    N_ext = sum(ext_sus_b)
    N_sus = nrow(storage_units)

    #  4.3 set variables
    #println("-- 4.3 set variables")
    @variables m begin
       (0 <=  SU_dispatch_fix[s=1:N_fix,t=1:T] <=
                (storage_units[fix_sus_b, :p_nom].*storage_units[fix_sus_b, :p_max_pu])[s])

        SU_dispatch_ext[s=1:N_ext,t=1:T] >= 0

        (0 <=  SU_store_fix[s=1:N_fix,t=1:T] <=
                 - (storage_units[fix_sus_b, :p_nom].*storage_units[fix_sus_b, :p_min_pu])[s])

        SU_store_ext[s=1:N_ext,t=1:T] >= 0

        SU_p_nom[s=1:N_ext] >= 0

        0 <= SU_soc_fix[s=1:N_fix,t=1:T] <= (storage_units[fix_sus_b,:max_hours]
                                            .*storage_units[fix_sus_b,:p_nom])[s]

        SU_soc_ext[s=1:N_ext,t=1:T] >= 0

        0 <=  SU_spill_fix[s=1:N_fix,t=1:T] <= inflow[:,fix_sus_b][t,s]

        0 <=  SU_spill_ext[s=1:N_ext,t=1:T] <= inflow[:,ext_sus_b][t,s]
    end

    # 4.4 set constraints for extendable storage_units
    #println("-- 4.4 set constraints for extendable storage_units")

    @constraints(m, begin
            su_dispatch_limit[s=1:N_ext,t=1:T], SU_dispatch_ext[s,t] <= SU_p_nom[s].*storage_units[ext_sus_b, :p_max_pu][s]
            su_storage_limit[s=1:N_ext,t=1:T], SU_store_ext[s,t] <= - SU_p_nom[s].*storage_units[ext_sus_b, :p_min_pu][s]
            su_capacity_limit[s=1:N_ext,t=1:T], SU_soc_ext[s,t] <= SU_p_nom[s].*storage_units[ext_sus_b, :max_hours][s]
    end)

    # 4.5 set charging constraint
    #println("-- 4.5 set charging constraint")
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
            su_logic1[s=is_cyclic_i], SU_soc[s,1] == (SU_soc[s,T]
                                        + storage_units[s,:efficiency_store] * SU_store[s,1]
                                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,1]
                                        + inflow[1,s] - SU_spill[s,1] )
            su_logic2[s=not_cyclic_i], SU_soc[s,1] == (storage_units[s,:state_of_charge_initial]
                                        + storage_units[s,:efficiency_store] * SU_store[s,1]
                                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,1]
                                        + inflow[1,s] - SU_spill[s,1])

            su_logic3[s=1:N_sus,t=2:T], SU_soc[s,t] == (SU_soc[s,t-1]
                                            + storage_units[s,:efficiency_store] * SU_store[s,t]
                                            - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,t]
                                            + inflow[t,s] - SU_spill[s,t] )

        end)

# --------------------------------------------------------------------------------------------------------

# 5. define stores
    println("Adding stores to the model.")

    # 5.1 set different stores types
    #println("-- 5.1 set different stores types")
    stores = network.stores

    fix_stores_b = .!stores[:e_nom_extendable]
    ext_stores_b = .!fix_stores_b

    inflow = get_switchable_as_dense(network, "stores", "inflow")

    # 5.2 iterator bounds
    #println("-- 5.2 iterator bounds")
    N_fix = sum(fix_stores_b)
    N_ext = sum(ext_stores_b)
    N_st = N_fix + N_ext

    #  5.3 set variables
    #println("-- 5.3 set variables")
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
    #println("-- 5.4 set constraints for extendable stores")

    @constraints(m, begin
            st_dispatch_limit[s=1:N_ext,t=1:T], ST_dispatch_ext[s,t] <= ST_e_nom[s].*stores[ext_stores_b, :e_max_pu][s]
            st_storage_limit[s=1:N_ext,t=1:T], ST_store_ext[s,t] <= - ST_e_nom[s].*stores[ext_stores_b, :e_min_pu][s]
            st_capacity_limit[s=1:N_ext,t=1:T], ST_soc_ext[s,t] <= ST_e_nom[s].*stores[ext_stores_b, :max_hours][s]
    end)

    # 5.5 set charging constraint
    #println("-- 5.5 set charging constraint")
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
            st_logic1[s=is_cyclic_i,t=1], ST_soc[s,t] == (ST_soc[s,T]
                                        + stores[s,:efficiency_store] * ST_store[s,t]
                                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]
                                        + inflow[t,s] - ST_spill[s,t])

            st_logic2[s=not_cyclic_i,t=1], ST_soc[s,t] == (stores[s,:state_of_charge_initial]
                                        + stores[s,:efficiency_store] * ST_store[s,t]
                                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]
                                        + inflow[t,s] - ST_spill[s,t])

            st_logic3[s=is_cyclic_i,t=2:T], ST_soc[s,t] == (ST_soc[s,t-1]
                                            + stores[s,:efficiency_store] * ST_store[s,t]
                                            - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]
                                            + inflow[t,s] - ST_spill[s,t])
        end)

# --------------------------------------------------------------------------------------------------------

# 6. power flow formulations

    println("Adding power flow formulation $formulation to the model.")

    # a.1 linear angles formulation
    if formulation == "angles_linear"

        # voltage angles
        @variable(m, THETA[1:N,1:T])

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

                == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) ))

        @constraint(m, flows[l=1:L, t=1:T], LN[l, t] == lines[:x_pu][l]^(-1) *     
                                                        ( THETA[busidx[lines[:bus0][l]], t]
                                                        - THETA[busidx[lines[:bus1][l]], t]
                                                        )
        )

        @constraint(m, slack[t=1:T], THETA[1,t] == 0 )

    elseif formulation == "angles_linear_integer_bigm"

        # voltage angles
        @variable(m, THETA[1:N,1:T])

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

                == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) ))

        for l=1:L
            for c in candidates[l]
                bigM_upper = maximum(candidates[l])*lines[:s_nom][l] + maximum(candidates[l])*lines[:x_pu][l]^(-1)*pi/6
                bigM_lower = maximum(candidates[l])*lines[:s_nom][l] + minimum(candidates[l])*lines[:x_pu][l]^(-1)*pi/6
                @constraints(m, begin 
                    [t=1:T],    (bigM_upper*(1-LN_opt[l,c]) 
                                            + (c*lines[:x_pu][l]^(-1) )
                                            *( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t]) - LN[l,t] >= 0)
                    [t=1:T],    (bigM_lower*(LN_opt[l,c]-1) 
                                            + (c*lines[:x_pu][l]^(-1) )
                                            *( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t]) - LN[l,t] <= 0)
                end)
            end
        end

        @constraint(m, slack[t=1:T], THETA[1,t] == 0 )

    # a.2 bilinear angles formulation (steps for x and s_nom separately specified)
    elseif formulation == "angles_bilinear_stepsspecs"

        # cannot be solved by Gurobi, needs Ipopt!
        # needs investment_type defined!

        # voltage angles
        @variable(m, THETA[1:N,1:T])

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

                == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) ))

        @NLconstraint(m, flows[l=1:L, t=1:T], LN[l, t] ==  (lines[:x_pu][l]^(-1) + lines[:x_step_pu][l]^(-1) * LN_inv[l]) *
                                                           (THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) )

        @constraint(m, slack[t=1:T], THETA[1,t] == 0)

    # a.3 bilinear angles formulation (steps derived from original s_nom and x)
    elseif formulation == "angles_bilinear"

        # cannot be solved by Gurobi, needs Ipopt!
        # needs investment_type defined!
        # requires starting reactance and line capacity!

        # voltage angles
        @variable(m, THETA[1:N,1:T])

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

                == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                   - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) ))

        @NLconstraint(m, flows[l=1:L, t=1:T], LN[l, t] ==  (1+LN_inv[l])*lines[:x_pu][l]^(-1) *
                                                           (THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) )

        @constraint(m, slack[t=1:T], THETA[1,t] == 0)

    # a.4 linear angles formulation (with McCormick relaxation)
    elseif formulation == "angles_relaxation_mccormick"

        # needs investment_type defined!

        # variables
        @variable(m, THETA[1:N,1:T])
        @variable(m, AUX_mu[1:L,1:T])
        @variable(m, AUX_xi[1:L,1:T])

        # load data in correct order
        loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

        @constraint(m, anglediff[l=1:L,t=1:T], AUX_xi[l,t] == THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t])

        max_angle_diff = pi / 6 
        min_angle_diff = -max_angle_diff

        @constraints(m, begin
            anglediff_upper_limit[l=1:L,t=1:T], AUX_xi[l,t] <= max_angle_diff
            anglediff_lower_limit[l=1:L,t=1:T], AUX_xi[l,t] >= min_angle_diff
        end )

        @constraint(m, voltages[n=1:N, t=1:T], 
        
                    sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                        .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                   - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) )

        @constraint(m, flows[l=1:L, t=1:T], 
            LN[l, t] == lines[:x_pu][l]^(-1) * AUX_xi[l,t]
                        + lines[:x_step_pu][l]^(-1) * AUX_mu[l,t]
        ) 

        inv_min = zeros(L)
        inv_max = zeros(L)
        for l=1:L
            inv_min[l] = (lines[:s_nom_min][l] - lines[:s_nom][l]) / lines[:s_nom_step][l]
            inv_max[l] = (lines[:s_nom_max][l] - lines[:s_nom][l]) / lines[:s_nom_step][l]
        end
        
        @constraints(m, begin
            mccormick1[l=1:L, t=1:T], (AUX_mu[l,t] >= inv_min[l] * AUX_xi[l,t]
                            + LN_inv[l] * min_angle_diff
                            - inv_min[l] * min_angle_diff)
            mccormick2[l=1:L, t=1:T], (AUX_mu[l,t] >= inv_max[l] * AUX_xi[l,t]
                            + LN_inv[l] * max_angle_diff
                            - inv_max[l] * max_angle_diff)
            mccormick3[l=1:L, t=1:T], (AUX_mu[l,t] <= inv_max[l] * AUX_xi[l,t]
                            + LN_inv[l] * min_angle_diff
                            - inv_max[l] * min_angle_diff)
            mccormick4[l=1:L, t=1:T], (AUX_mu[l,t] <= inv_min[l] * AUX_xi[l,t]
                            + LN_inv[l] * max_angle_diff
                            - inv_min[l] * max_angle_diff)
        end)

        @constraint(m, slack[t=1:T], THETA[1,t] == 0)

    # a.5 linear angles formulation (with relaxation following Taylor2015a)
    elseif formulation == "angles_relaxation_taylor2015a"

        # not working yet

        inv_max = zeros(L)
        for l=1:L
            inv_max[l] = (lines[:s_nom_max][l] - lines[:s_nom][l]) / lines[:s_nom_step][l]
        end

        # variables
        @variable(m, THETA[1:N,1:T])
        @variable(m, AUX_XI[1:L,1:T])

        # load data in correct order
        loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

        # (2)
        @constraint(m, voltages[n=1:N, t=1:T], (

                    sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                        .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                   - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) ))

        # (2)
        @constraint(m, flows[l=1:L, t=1:T], LN[l, t] == lines[:x_pu][l]^(-1) *
                                                        ( THETA[busidx[lines[:bus0][l]], t]
                                                        - THETA[busidx[lines[:bus1][l]], t]
                                                        )

                                                      + AUX_XI[l,t] )

        @constraint(m, slack[t=1:T], THETA[1,t] == 0)

        # (1)
        L_existing = Int64[]
        for l=1:L
            if lines[:s_nom][l] > 0
                push!(L_existing, l)
            end
        end

        @constraints(m, begin
            [l in L_existing, t=1:T], lines[:x_pu][l]^(-1) *
                            ( THETA[busidx[lines[:bus0][l]], t]
                            - THETA[busidx[lines[:bus1][l]], t]
                            ) <= lines[:s_nom][l]

            [l in L_existing, t=1:T], lines[:x_pu][l]^(-1) *
                            ( THETA[busidx[lines[:bus1][l]], t]
                            - THETA[busidx[lines[:bus0][l]], t]
                            ) >= lines[:s_nom][l]
            end
        )

        # (4)
        @constraints(m, begin
            [l=1:L, t=T],  AUX_XI[l,t] <= lines[:s_nom_step][l] * LN_inv[l]
            [l=1:L, t=T],  AUX_XI[l,t] >= (-1) * lines[:s_nom_step][l] * LN_inv[l]
        end )

        # calculate Ms
        line_paths = get_line_paths(network)
        m_factor = zeros(L)
        for l=1:L
            for k in line_paths[l]
                m_factor[l] += lines[:s_nom_step][k] * lines[:x_step_pu][k]
            end
        end

        # (3)
        @constraints(m, begin
            [l=1:L, t=T],  AUX_XI[l,t] <= lines[:x_step_pu][l]^(-1) * LN_inv[l] * m_factor[l]
            [l=1:L, t=T],  AUX_XI[l,t] >= (-1) * lines[:x_step_pu][l]^(-1) * LN_inv[l] * m_factor[l]
        end )

        # (5)
        @constraints(m, begin
            [l=1:L, t=T],  AUX_XI[l,t] - lines[:x_step_pu][l]^(-1) * inv_max[l] * ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t]) <= lines[:x_step_pu][l]^(-1) * m_factor[l] * (inv_max[l] - LN_inv[l])
            [l=1:L, t=T],  AUX_XI[l,t] - lines[:x_step_pu][l]^(-1) * inv_max[l] * ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t]) >= (-1) * ( lines[:x_step_pu][l]^(-1) * m_factor[l] * (inv_max[l] - LN_inv[l]) )
        end )

    # b.1 linear kirchhoff formulation
    elseif formulation == "kirchhoff_linear"

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

    # b.2 bilinear kirchhoff formulation (steps derived from original s_nom and x)
    elseif formulation == "kirchhoff_bilinear"

        # cannot be solved by Gurobi, needs Ipopt!
        # needs investment_type defined!
        # requires starting reactance and line capacity!

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
                @NLconstraint(m, line_cycle_constraint[c=1:length(cycles_branch), t=1:T] ,
                        sum(      directions[c][l] 
                                * lines[cycles_branch[c], :x_pu][l] 
                                * (1+LN_inv[cycles_branch[c]][l])^(-1)
                                * LN[cycles_branch[c],t][l]
                           for l=1:length(directions[c])
                           )               
                            == 0)
            # elseif attribute==:r
            #     @constraint(m, link_cycle_constraint[c=1:length(cycles_branch), t=1:T] ,
            #             dot(directions[c] .* links[cycles_branch[c], :r]/380. , LK[cycles_branch[c],t]) == 0)
            end
        end

    # c.1 linear ptdf formulation
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

# 7. set global_constraints

    println("Adding global CO2 constraints to the model.")

    if nrow(network.global_constraints)>0 && in("primary_energy", network.global_constraints[:type])
        co2_limit = network.global_constraints[network.global_constraints[:name].=="co2_limit", :constant]
        nonnull_carriers = network.carriers[network.carriers[:co2_emissions].!=0, :]
        carrier_index(carrier) = findin(generators[:carrier], [carrier])
        @constraint(m, sum(sum(dot(1./generators[carrier_index(carrier) , :efficiency],
                    G[carrier_index(carrier),t]) for t=1:T)
                    * select_names(network.carriers, [carrier])[:co2_emissions]
                    for carrier in network.carriers[:name]) .<=  co2_limit)
    end

# --------------------------------------------------------------------------------------------------------

# 8. set objective function
    println("Adding objective to the model.")

    @objective(m, Min,
                          sum(dot(generators[:marginal_cost], G[:,t]) for t=1:T)
                        + dot(generators[ext_gens_b,:capital_cost], G_p_nom[:] )
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

    println("Finished building model.")

    return m

end
