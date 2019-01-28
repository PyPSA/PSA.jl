using JuMP
using Test


include("auxilliaries.jl")


function lopf_pathway(n, solver; extra_functionality=nothing, 
    investment_period=nothing, invest_at_first_sn=false, fix_p_nom=nothing)
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


    m = Model(solver=solver)  # from JuMP

    calculate_dependent_values!(n)   # function in auxilliaries
    buses = n.buses.axes[1].val
    busidx = idx(n.buses)    # Dict{Busname, Index}
    reverse_busidx = rev_idx(n.buses)
    N, = size(n.buses)
    T, = size(n.snapshots) #normally snapshots
    #t_ip = index of snapshot, when you invest, IP number of investments
    t_ip, IP = get_investment_periods(n, investment_period, invest_at_first_sn)
    @info("Solve network with $IP investment timesteps.")
    total_capacity_line_ext = 1e7
    @info("The total capacity of the extendable lines in all investment periods is $total_capacity_line_ext .")
    (size(n.loads_t.p)[1] != T ?
        replace_attribute!(n, :loads_t, :p, n.loads_t.p_set) : nothing)

# --------------------------------------------------------------------------------------------------------

# 1. add all generators to the model
    # 1.1 set different generator types
    generators = n.generators

    ext_gens_b = BitArray(n.generators[:, "p_nom_extendable"])
    com_gens_b = BitArray(n.generators[:, "commitable"])
    fix_gens_b = (ext_gens_b .| com_gens_b); fix_gens_b = .!fix_gens_b #two steps is faster

    # 1.2 fix iterating and variable bounds
    N_fix = sum(fix_gens_b)
    N_ext = sum(ext_gens_b)
    N_com = sum(com_gens_b)

    p_max_pu = get_switchable_as_dense(n, "generators", "p_max_pu") # matrix snapshot x generator["p_max_pu"]
    p_min_pu = get_switchable_as_dense(n, "generators", "p_min_pu")

    @test all(broadcast(<, p_min_pu, p_max_pu)) == true

    p_nom = float.(n.generators[:, "p_nom"]) # Array{Any,1} formatted to Array{Float64,1}

    Ub_fix = p_max_pu[:,fix_gens_b] .* p_nom[fix_gens_b]'
    Lb_fix = p_min_pu[:,fix_gens_b] .* p_nom[fix_gens_b]'
    
    Ub_com = p_max_pu[:,com_gens_b] .* p_nom[com_gens_b]'
    Lb_com = p_min_pu[:,com_gens_b] .* p_nom[com_gens_b]'

    Ub_ext_nom = n.generators[ext_gens_b, "p_nom_max"] # Array(number generators x 1)  p_nom_max
    Lb_ext_nom = n.generators[ext_gens_b, "p_nom_min"]

    # 1.3 add generator variables to the model
    @variables m begin

        Lb_fix[t,gr] <= G_fix[t=1:T, gr=1:N_fix] <= Ub_fix[t,gr] 

                        G_ext[t=1:T, gr=1:N_ext]
        
                        G_status[t=1:T, gr=1:N_com]

        0            <= G_p_nom_ext[ip=1:IP, gr=1:N_ext] <= Ub_ext_nom[gr]

        0            <= G_p_nom_red[ip=1:IP, gr=1:N_ext] <= Ub_ext_nom[gr]

        Lb_com[t,gr] <= G_com[t=1:T,gr=1:N_com] <= Ub_com[t,gr]
    end

    # 1.4 add time-varying upper bound for p_nom
    @expression(m, G_p_nom[t=1:T, gr=1:N_ext], 0)

    for gr=1:N_ext
        @assert p_nom[ext_gens_b][gr] <= Ub_ext_nom[gr]
        @assert p_nom[ext_gens_b][gr] >= Lb_ext_nom[gr]
    end
    aggregate_investments!(G_p_nom, p_nom[ext_gens_b], 
                           G_p_nom_ext - G_p_nom_red, t_ip, true)

    @constraints(m, begin
        [t=1:T, gr=1:N_ext], G_p_nom[t, gr]  >= Lb_ext_nom[gr]
        [t=1:T, gr=1:N_ext], G_p_nom[t, gr]  <= Ub_ext_nom[gr]
    end)


    # # add fixed p_nom at a certain snapshot
    if length(n.generators_t.p_nom_opt)!=0
        @info("There are fixed p_nom values for the generators.")
        p_nom_set = .!isnan.(n.generators_t.p_nom_opt)
        val = n.generators_t.p_nom_opt[p_nom_set]
        index = findn(p_nom_set)
        for (t,gr, v) in zip(index...,val)
            @constraint(m,G_p_nom[t, gr] == v)
        end

    end

    # 1.5 set constraints for generator

    Ub_ext = p_max_pu[:,ext_gens_b] .* G_p_nom
    Lb_ext = p_min_pu[:,ext_gens_b] .* G_p_nom  

    @constraints(m, begin
        [t=1:T,gr=1:N_ext], G_ext[t,gr] >= Lb_ext[t,gr]
        [t=1:T,gr=1:N_ext], G_ext[t,gr] <= Ub_ext[t,gr]

        [t=1:T,gr=1:N_com], G_com[t,gr] - G_com[t,gr].*G_status[t,gr] == 0
    end)


    G = [G_fix G_ext G_com] # G is the concatenated variable array
    # sort generators the same way
    generators = generators[[findall(fix_gens_b) ; findall(ext_gens_b); findall(com_gens_b)],:] 
    # new booleans
    ext_gens_b = BitArray(generators[:, "p_nom_extendable"])
    com_gens_b = BitArray(generators[:, "commitable"])
    fix_gens_b = (ext_gens_b .| com_gens_b); fix_gens_b = .!fix_gens_b 

    p_max_pu = 0; p_min_pu = 0

# --------------------------------------------------------------------------------------------------------

# 2. add all lines to the model
    # 2.1 set different lines types
    lines = n.lines
    ext_lines_b = BitArray(lines[:, "s_nom_extendable"])
    fix_lines_b = .!ext_lines_b
    s_nom_min_lines_b = BitArray(lines[:, "s_nom_min"] .> 0)

    # 2.2 iterator bounds
    N_fix = sum(fix_lines_b)
    N_ext = sum(ext_lines_b)

    Lb_fix = -lines[fix_lines_b,"s_nom"]
    Ub_fix = lines[fix_lines_b,"s_nom"]
    
    Lb_ext_nom = lines[ext_lines_b,"s_nom_min"]
    Ub_ext_nom = lines[ext_lines_b,"s_nom_max"]

    s_nom = float.(n.lines[:, "s_nom"])

    # 2.3 add variables
    @variables m begin
        Lb_fix[l]  <=  LN_fix[t=1:T,l=1:N_fix] <= Ub_fix[l]

        LN_ext[t=1:T,l=1:N_ext]
        
        LN_s_nom_ext[ip=1:IP, l=1:N_ext] >=  0
    end


    # 2.4 add time-varying upper bound for s_nom
    @expression(m, LN_s_nom[t=1:T, l=1:N_ext], 0)

    for l=1:N_ext
        @assert s_nom[ext_lines_b][l] <= Ub_ext_nom[l]
        @assert s_nom[ext_lines_b][l]  >= Lb_ext_nom[l]
    end

    aggregate_investments!(LN_s_nom, s_nom[ext_lines_b], LN_s_nom_ext, t_ip, false)

    # 2.4 add line constraint for extendable lines
    @constraints(m, begin
            [t=1:T,l=1:N_ext], LN_s_nom[t,l] <=  Ub_ext_nom[l]
            [t=1:T,l=1:N_ext], LN_s_nom[t,l] >= Lb_ext_nom[l]
            [ip=1:IP], sum(LN_s_nom_ext[ip, :]) <= total_capacity_line_ext/IP 
    end)

    Ub_ext = Ub_ext_nom' .* LN_s_nom

    @constraints(m, begin
        [t=1:T, l=1:N_ext], LN_ext[t,l] >= -LN_s_nom[t,l]
        [t=1:T, l=1:N_ext], LN_ext[t,l] <= LN_s_nom[t,l]
    end )

    LN = [LN_fix LN_ext]
    lines = lines[[findall(fix_lines_b) ; findall(ext_lines_b)], : ]

    ext_lines_b = BitArray(lines[:, "s_nom_extendable"])
    fix_lines_b = .! ext_lines_b
    @info("lines passed")
# ------------------------------------------------------------------------------

# 3. add all links to the model
    # 3.1 set different link types
    links = n.links
    ext_links_b = BitArray(links[:,"p_nom_extendable"])
    fix_links_b = .!ext_links_b

    # 3.2 iterator bounds
    N_fix = sum(fix_links_b)
    N_ext = sum(ext_links_b)

    Lb_fix = links[fix_links_b, "p_nom"] .* links[fix_links_b, "p_min_pu"]
    Ub_fix = links[fix_links_b, "p_nom"] .* links[fix_links_b, "p_max_pu"]
    
    Lb_nom = links[ext_links_b, "p_nom_min"]
    Ub_nom = links[ext_links_b, "p_nom_max"]
    
    p_nom = float.(n.links[:, "p_nom"])
    
    #  3.3 set link variables
    @variables m begin
        Lb_fix[l]  <=  LK_fix[t=1:T,l=1:N_fix]  <= Ub_fix[l]

        LK_ext[t=1:T,l=1:N_ext]

        LK_p_nom_ext[ip=1:IP, l=1:N_ext] >= 0
    end

    @expression(m, LK_p_nom[t=1:T, l=1:N_ext], 0)
    for l=1:N_ext
        @assert p_nom[ext_links_b][l] <= Ub_ext_nom[l]
        @assert p_nom[ext_links_b][l]  >= Lb_ext_nom[l]
    end
    aggregate_investments!(LK_p_nom, p_nom[ext_links_b], LK_p_nom_ext, t_ip, true)


    # 3.4 set constraints for extendable links
    @constraints(m, begin
    [t=1:T,l=1:N_ext], LK_ext[t,l] >= Lb_nom[l]
    [t=1:T,l=1:N_ext], LK_ext[t,l] <= LK_p_nom[t,l]
    end)

    @constraints(m, begin
    [t=1:T,l=1:N_ext], LK_p_nom[t,l] >= Lb_nom[l]
    [t=1:T,l=1:N_ext], LK_p_nom[t,l] <= Ub_nom[l]
    [ip=1:IP], sum(LK_p_nom_ext[ip, :]) <= total_capacity_line_ext/IP 
    # change later with cost limit for DC lines
    end)

    LK = [LK_fix LK_ext]
    links = links[[findall(fix_links_b) ; findall(ext_links_b)], : ]
    ext_links_b = BitArray(links[:, "p_nom_extendable"])
    fix_links_b = .!ext_links_b

    @info("links passed")
# ------------------------------------------------------------------------------
# 4. define storage_units
    # 4.1 set different storage_units types
    storage_units = n.storage_units
    fix_sus_b = .!BitArray(storage_units[:, "p_nom_extendable"])
    ext_sus_b = BitArray(.!fix_sus_b)

    inflow = get_switchable_as_dense(n, "storage_units", "inflow")

    # 4.2 iterator bounds
    N_fix = sum(fix_sus_b)
    N_ext = sum(ext_sus_b)
    N_sus = N_fix + N_ext
    
    Ub_dist_fix = storage_units[fix_sus_b, "p_nom"].*storage_units[fix_sus_b, "p_max_pu"]
    Ub_stor_fix = storage_units[fix_sus_b, "p_nom"].*storage_units[fix_sus_b, "p_min_pu"] 
    Ub_soc_fix = storage_units[fix_sus_b, "max_hours"] .* storage_units[fix_sus_b, "p_nom"]
    Ub_spill_fix = inflow[:,fix_sus_b]
    Ub_spill_ext = inflow[:,ext_sus_b]

    #  4.3 set variables
    @variables m begin
        0 <=  SU_dispatch_fix[t=1:T,s=1:N_fix] <= Ub_dist_fix[s]

        SU_dispatch_ext[t=1:T,s=1:N_ext] >= 0

        0 <=  SU_store_fix[t=1:T,s=1:N_fix] <= - Ub_stor_fix[s]
        
        SU_store_ext[t=1:T,s=1:N_ext] >= 0

        SU_p_nom[s=1:N_ext] >= 0

        0 <= SU_soc_fix[t=1:T,s=1:N_fix] <= Ub_soc_fix[s]
        
        SU_soc_ext[t=1:T,s=1:N_ext] >= 0

        0 <=  SU_spill_fix[t=1:T,s=1:N_fix] <= Ub_spill_fix[t,s]
        
        0 <=  SU_spill_ext[t=1:T,s=1:N_ext] <= Ub_spill_ext[t,s]
        end

    # 4.4 set constraints for extendable storage_units
    Ub_dist_ext = SU_p_nom .* storage_units[ext_sus_b, "p_max_pu"]
    Ub_stor_ext = - SU_p_nom .* storage_units[ext_sus_b, "p_min_pu"]
    Ub_soc_ext = SU_p_nom .* storage_units[ext_sus_b, "max_hours"]

    @constraints(m, begin
            [t=1:T,s=1:N_ext], SU_dispatch_ext[t,s] <= Ub_dist_ext[s]
            [t=1:T,s=1:N_ext], SU_store_ext[t,s] <= Ub_stor_ext[s]
            [t=1:T,s=1:N_ext], SU_soc_ext[t,s] <= Ub_soc_ext[s]
    end)
    
    # 4.5 set charging constraint
    SU_dispatch = [SU_dispatch_fix SU_dispatch_ext]
    SU_store = [SU_store_fix SU_store_ext]
    SU_soc = [SU_soc_fix SU_soc_ext]
    SU_spill = [SU_spill_fix SU_spill_ext]

    storage_units = storage_units[[findall(fix_sus_b) ; findall(ext_sus_b)], :]
    inflow = inflow[:,[findall(fix_sus_b); findall(ext_sus_b)]]
    ext_sus_b = BitArray(storage_units[:,"p_nom_extendable"])

    
    is_cyclic_i = findall(storage_units[:, "cyclic_state_of_charge"])
    not_cyclic_i = findall(.!BitArray(storage_units[:,"cyclic_state_of_charge"]))
  
    
    SU_soc_change = (storage_units[:,"efficiency_store"]' .* SU_store 
                    - 1 ./storage_units[:,"efficiency_dispatch"]' .* SU_dispatch 
                    + inflow - SU_spill)


    @constraints(m, begin
            [s=is_cyclic_i], SU_soc[1,s] == SU_soc[T,s] + SU_soc_change[1,s]
            [s=not_cyclic_i], SU_soc[1,s] == storage_units[s,"state_of_charge_initial"] + SU_soc_change[1,s] 
            [t=2:T,s=1:N_sus], SU_soc[t,s] == SU_soc[t-1,s] + SU_soc_change[t,s]
        end)
    @info("storage units passed")
# ------------------------------------------------------------------------------

# 5. define stores
# 5.1 set different stores types
    stores = n.stores
    fix_stores_b = .!BitArray(stores[:,"e_nom_extendable"])
    ext_stores_b = BitArray(.!fix_stores_b)

    inflow = get_switchable_as_dense(n, "stores", "inflow")

    # 5.2 iterator bounds
    N_fix = length(fix_stores_b) == 0 ? 0 : sum(fix_stores_b)
    N_ext = length(ext_stores_b) == 0 ? 0 : sum(ext_stores_b)
    N_st = N_fix + N_ext

    Ub_dist_fix = stores[fix_stores_b, "e_nom"].*stores[fix_stores_b, "e_max_pu"]
    Ub_stor_fix = stores[fix_stores_b, "e_nom"].*stores[fix_stores_b, "e_min_pu"] 
    Ub_soc_fix = stores[fix_stores_b, "max_hours"] .* stores[fix_stores_b, "e_nom"]
    Ub_spill_fix = inflow[:,fix_stores_b]
    Ub_spill_ext = inflow[:,ext_stores_b]

    #  5.3 set variables
    @variables m begin
        0 <=  ST_dispatch_fix[t=1:T,s=1:N_fix] <= Ub_dist_fix[s]

        ST_dispatch_ext[t=1:T,s=1:N_ext] >= 0

        0 <=  ST_store_fix[t=1:T,s=1:N_fix] <= - Ub_stor_fix[s]
        
        ST_store_ext[t=1:T,s=1:N_ext] >= 0

        ST_e_nom[s=1:N_ext] >= 0

        0 <= ST_soc_fix[t=1:T,s=1:N_fix] <= Ub_soc_fix[s]
        
        ST_soc_ext[t=1:T,s=1:N_ext] >= 0

        0 <=  ST_spill_fix[t=1:T,s=1:N_fix] <= Ub_spill_fix[t,s]
        
        0 <=  ST_spill_ext[t=1:T,s=1:N_ext] <= Ub_spill_ext[t,s]
        end


    # 5.4 set constraints for extendable stores
    Ub_dist_ext = ST_e_nom .* stores[ext_stores_b, "e_max_pu"]
    Ub_stor_ext = - ST_e_nom .* stores[ext_stores_b, "e_min_pu"]
    Ub_soc_ext = ST_e_nom .* stores[ext_stores_b, "max_hours"]

    @constraints(m, begin
            [t=1:T,s=1:N_ext], ST_dispatch_ext[t,s] <= Ub_dist_ext[s]
            [t=1:T,s=1:N_ext], ST_store_ext[t,s] <= Ub_stor_ext[s]
            [t=1:T,s=1:N_ext], ST_soc_ext[t,s] <= Ub_soc_ext[s]
    end)

    # 5.5 set charging constraint
    ST_dispatch = [ST_dispatch_fix ST_dispatch_ext]
    ST_store = [ST_store_fix ST_store_ext]
    ST_soc = [ST_soc_fix ST_soc_ext]
    ST_spill = [ST_spill_fix ST_spill_ext]

    stores = stores[[findall(fix_stores_b) ; findall(ext_stores_b)], :]
    inflow = inflow[:,[findall(fix_stores_b); findall(ext_stores_b)]]
    ext_stores_b = BitArray(stores[:,"e_nom_extendable"])


    is_cyclic_i = findall(stores[:, "cyclic_state_of_charge"])
    not_cyclic_i = findall(.!stores[:,"cyclic_state_of_charge"])

    ST_soc_change = (stores[:,"efficiency_store"]' .* ST_store 
                    - 1 ./stores[:,"efficiency_dispatch"]' .* ST_dispatch 
                    + inflow - ST_spill)


    @constraints(m, begin
            [s=is_cyclic_i], ST_soc[1,s] == ST_soc[T,s] + ST_soc_change[1,s]
            [s=not_cyclic_i], ST_soc[1,s] == stores[s,"state_of_charge_initial"] + ST_soc_change[1,s] 
            [t=2:T,s=1:N_st], ST_soc[t,s] == ST_soc[t-1,s] + ST_soc_change[t,s]
        end)

    tic()

    @info("stores passed")
# ------------------------------------------------------------------------------

## 6. define nodal balance constraint

    #load data in correct order
    loads = n.loads_t.p.data

    add_missing_buses(d) = (for b=setdiff(buses, keys(d)) d[b]  = Int[] end; d)

    g_bus = string.(generators[:, "bus"]) |> grouped_array |> add_missing_buses
    l_bus = string.(n.loads[:, "bus"]) |> grouped_array |> add_missing_buses
    ln_bus0 = string.(lines[:, "bus0"]) |> grouped_array |> add_missing_buses
    ln_bus1 = string.(lines[:, "bus1"]) |> grouped_array |> add_missing_buses
    lk_bus0 = string.(links[:, "bus0"]) |> grouped_array |> add_missing_buses
    lk_bus1 = string.(links[:, "bus1"]) |> grouped_array |> add_missing_buses
    su_bus = string.(storage_units[:, "bus"]) |> grouped_array |> add_missing_buses

    linkefficieny = float.(links[:, "efficiency"])
    
    @constraint(m, balance[b=1:N, t=1:T], (
          sum(G[t, g_bus[buses[b]]])
        + sum(LN[t, ln_bus1[buses[b]] ])
        + zsum(linkefficieny[lk_bus1[buses[b]]] .* LK[t, lk_bus1[buses[b]]] )
        + sum(SU_dispatch[t, su_bus[buses[b]]])

        - sum(loads[t, l_bus[buses[b]]])
        - sum(LN[t, ln_bus0[buses[b]] ])
        - sum(LK[t, lk_bus0[buses[b]] ])
        - sum(SU_store[t, su_bus[buses[b]]])

          == 0 ))

    @info("nodal balance passed")
# --------------------------------------------------------------------------------------------------------

# 7. set Kirchhoff Voltage Law constraint

# since cyclebasis is not yet supported in LightGraphs, use the pyimported
# netwrokx in order to define all cycles. The cycle_basis returns a list of
# cycles, indicating the connected buses. For each cycle the connecting branches
# and their directions (signs) have to be determined.

# Might be nessecary to loop over all subgraphs as
# for (sn, sub) in enumerate(weakly_connected_components(g))
#     g_sub = induced_subgraph(g, sub)[1]

    # @variable(m, y, Bin)   # add binary variable
    # lines_used = LN .!= 0
    (branches, var, attribute) = (lines[fix_lines_b, :], LN_fix, "x_pu")
    @show(size(branches))
    @info("Kirchhoff 1 passed")
    # Here is the problematic function
    cycles, dirs = get_directed_cycles(n, branches)
    @info("Kirchhoff 2 passed")
    @constraint(m, line_cycle_constraint[c=1:length(cycles), t=1:T] ,
            dot(dirs[c] .* float.(branches[cycles[c], attribute]) .* 1e5 ,
                var[t, cycles[c]]) == 0)
    @info("Kirchhoff 3 passed")
# --------------------------------------------------------------------------------------------------------

# 8. set global_constraints
# only for co2_emissions till now

    if size(n.global_constraints)[1]>0# && in("primary_energy", n.global_constraints[:,"type"])
        co2_limit = n.global_constraints["CO2Limit", "constant"]
        nonnull_carriers = n.carriers[n.carriers[:,"co2_emissions"].!=0, :]

        carrier_index(carrier) = findall(in(string.(generators[:,"carrier"])), [carrier])

        @constraint(m, sum(sum(dot(1/float.(generators[carrier_index(carrier) , "efficiency"]),
                    G[t,carrier_index(carrier)]) for t=1:T)
                    * n.carriers[carrier, "co2_emissions"]
                    for carrier in nonnull_carriers.axes[1].val) <=  co2_limit)
    end

    @info("global constraints passed")

# --------------------------------------------------------------------------------------------------------
# muss noch umgeschrieben werden, maintenance_cost für alle Generatoren, nicht nur extenable
# maintenance_cost = fix costs per hour
# marginal_cost = variable cost per hour
# 10. set objective function

    @objective(m, Min,

            zsum(float.(generators[:,"marginal_cost"])'.* G .*n.snapshot_weightings)
            
            + zsum(float.(generators[ext_gens_b,"maintenance_cost"])'.* G_p_nom .*n.snapshot_weightings)

            + zsum(float.(generators[ext_gens_b,"capital_cost"])' .* G_p_nom_ext)
            + zsum(float.(generators[ext_gens_b,"capital_cost"])' .*0.01 .* G_p_nom_red)

            + zsum(float.(lines[ext_lines_b,"capital_cost"])' .* LN_s_nom)
            + zsum(float.(links[ext_links_b,:"capital_cost"])' .*  LK_p_nom) 

            + zsum(float.(storage_units[:,"marginal_cost"])' .* SU_dispatch .*n.snapshot_weightings) 
            + zdot(float.(storage_units[ext_sus_b, "capital_cost"]), SU_p_nom)

            + zsum(float.(stores[:,"marginal_cost"])' .* ST_dispatch .*n.snapshot_weightings) 
            + zdot(float.(stores[ext_stores_b, "capital_cost"]), ST_e_nom[:])
        )

    status = solve(m)

# --------------------------------------------------------------------------------------------------------

# 11. extract optimisation results
    if status==:Optimal
        order_components = get_comp_order(n)  # get the right index order
        # generators
        extract_opt_result!(n, G_p_nom, :generators_t, "p_nom_opt", n.generators.axes[1].val[ext_gens_b])
        extract_opt_result!(n, G, :generators_t, "p", n.generators.axes[1].val)
        n.generators = reindex(n.generators, index = n.generators.axes[1].val)
        @info("Passed")

        # lines
        orig_line_order = n.lines.axes[1].val
        n.lines[fix_lines_b,"s_nom_opt"] = n.lines[fix_lines_b, "s_nom"]
        extract_opt_result!(n, LN, :lines_t, "p0", n.lines.axes[1].val)
        extract_opt_result!(n, LN_s_nom, :lines_t, "s_nom_opt", n.lines.axes[1].val[ext_lines_b])
        n.lines = reindex(n.lines, index = orig_line_order)

        # links
        orig_link_order = n.links.axes[1].val
        n.links[fix_links_b,"p_nom_opt"] = n.links[fix_links_b, "p_nom"]
        extract_opt_result!(n, LK, :links_t, "p0", n.links.axes[1].val)
        extract_opt_result!(n, LK_p_nom, :links_t, "p_nom_opt", n.links.axes[1].val[ext_links_b])
        n.links = reindex(n.links, index = orig_link_order)

        # storage_units
        orig_sus_order = n.storage_units.axes[1].val
        n.storage_units = storage_units
        n.storage_units[fix_sus_b,"p_nom_opt"] = n.storage_units[fix_sus_b, "p_nom"]
        n.storage_units[ext_sus_b,"p_nom_opt"] = getvalue(SU_p_nom)
        extract_opt_result!(n, SU_spill, :storage_units_t, "spill", n.storage_units.axes[1].val)
        data = AxisArray(getvalue(SU_dispatch .- SU_store), 
                                                Axis{:row}(n.snapshots), 
                                                Axis{:col}(n.storage_units.axes[1].val))
        replace_attribute!(n, :storage_units_t, "p", data)
        data = AxisArray(getvalue(SU_soc), Axis{:row}(n.snapshots), 
                                           Axis{:col}(n.storage_units.axes[1].val))
        replace_attribute!(n, :storage_units_t, "state_of_charge", data)
        n.storage_units = reindex(n.storage_units, index = orig_sus_order)
        @info("Passed")

        # stores
        orig_stores_order = n.stores.axes[1].val
        n.stores = stores
        n.stores[fix_stores_b,"e_nom_opt"] = n.stores[fix_stores_b, "e_nom"]
        n.stores[ext_stores_b,"e_nom_opt"] = getvalue(ST_e_nom)
        data = AxisArray(getvalue(ST_spill), Axis{:row}(n.snapshots), 
                                             Axis{:col}(n.stores.axes[1].val))
        replace_attribute!(n, :stores_t, "spill", data)
        data = AxisArray(getvalue(getvalue(ST_dispatch .- ST_store)), 
                                                Axis{:row}(n.snapshots), 
                                                Axis{:col}(n.stores.axes[1].val))
        replace_attribute!(n, :stores_t, "e", data)                   
        data = AxisArray(getvalue(ST_soc), Axis{:row}(n.snapshots), 
                                           Axis{:col}(n.stores.axes[1].val))
        replace_attribute!(n, :stores_t, "state_of_charge", data)
        n.stores = reindex(n.stores, index = orig_stores_order)

        
        # buses TODO : need reindexing along nonexisting columns
        # n.buses_t["p"] =  group(n.lines_t["p0"], n.lines[:,"bus0"], sum,  axis=2) + 
        #                   group(n.links_t["p0"], n.links[:,"bus0"], sum,  axis=2) - 
        #                   group(n.lines_t["p0"], n.lines[:,"bus1"], sum,  axis=2) -
        #                   group(n.links_t["p0"], n.links[:,"bus1"], sum,  axis=2)
        
        align_component_order!(n)
        println("Reduce cost to $(m.objVal)")
    end
    return m
end
