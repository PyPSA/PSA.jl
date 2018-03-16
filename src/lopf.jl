import Gurobi
using JuMP
# using Ipopt


include("auxilliaries.jl")


function lopf(network; solver_options...)
    env = Gurobi.Env()
    Gurobi.setparams!(env; solver_options...)
    m = Model(solver=Gurobi.GurobiSolver())
    
    calculate_dependent_values!(network)
    buses = network.buses
    reverse_busidx = rev_idx(buses)
    busidx = idx(buses)
    N = nrow(network.buses)
    T = nrow(network.snapshots) #normally snapshots
    nrow(network.loads_t["p"])!=T ? network.loads_t["p"]=network.loads_t["p_set"] : nothing


# 1. add all generators to the model
    # 1.1 set different generator types
    fix_gens_b = ((.!network.generators[:p_nom_extendable]) .& (.!network.generators[:commitable]))
    ext_gens_b = convert(BitArray, network.generators[:p_nom_extendable])
    com_gens_b = convert(BitArray, network.generators[:commitable])
    generators = network.generators

    # 1.2 fix bounds for iterating
    G_fix = sum(fix_gens_b)
    G_ext = sum(ext_gens_b)
    G_com = sum(com_gens_b)

    p_max_pu = get_switchable_as_dense(network, "generators", "p_max_pu")
    p_min_pu = get_switchable_as_dense(network, "generators", "p_min_pu")

   # 1.3 add generator variables to the model
    @variables m begin

        (generators[fix_gens_b,:p_nom][gr]*p_min_pu[:,fix_gens_b][t,gr] <= g_fix[gr=1:G_fix,t=1:T]
                                <= generators[fix_gens_b,:p_nom][gr]*p_max_pu[:,fix_gens_b][t,gr])

        g_ext[gr=1:G_ext,t = 1:T]
        generators[ext_gens_b, :p_nom_min][gr] <=  gen_p_nom[gr=1:G_ext] <= generators[ext_gens_b,:p_nom_max][gr]

        g_status[gr=1:G_com,t=1:T]
        (generators[com_gens_b,:p_nom][gr].*p_min_pu[:,com_gens_b][t,gr] <= g_com[gr=1:G_com,t=1:T]
                                <= generators[com_gens_b,:p_nom][gr]*p_max_pu[:,com_gens_b][t,gr])
    end

    # 1.4 set constraints for generators

    @constraints(m, begin
        [gr=1:G_ext,t=1:T], g_ext[gr,t] >= gen_p_nom[gr]*p_min_pu[:,ext_gens_b][t,gr]
        [gr=1:G_ext,t=1:T], g_ext[gr,t] <= gen_p_nom[gr]*p_max_pu[:,ext_gens_b][t,gr]

        [gr=1:G_com,t=1:T], g_com[gr,t] - g_com[gr,t].*g_status[gr,t] == 0
    end)


    gn = [g_fix; g_ext; g_com] # gn is the concatenated variable array
    generators = [generators[fix_gens_b,:]; generators[ext_gens_b,:]; generators[com_gens_b,:] ] # sort generators the same

# commitable still to work on:
    # append_idx_col!(generators)
    # g_up_time_b = generators[(com_gens_b .& generators[:min_up_time].>0), :idx]
    # g_down_time_b = generators[(com_gens_b .& generators[:min_down_time].>0), :idx]

    # @constraints(m, begin
        # [gr=G_com,t=1], (sum(g_status[gr,j] for j=t:min.(t+generators_com[gr,:min_up_time]-1,T))
        #                     >=
        #                     # generators_com[gr,:min_up_time].*g_status[gr,t]
        #                     generators_com[gr,:min_up_time].*generators_com[gr, :initial_status])
        # [gr=G_com,t=2:T], (sum(g_status[gr,j] for j=t:min.(t+generators_com[gr,:min_up_time]-1,T))
        #                     >= generators_com[gr,:min_up_time].*g_status[gr,t]
        #                     - generators_com[gr,:min_up_time].*g_status[gr,t-1])
        #
        # [gr=G_com,t=1], (generators_com[gr,:min_down_time]
        #                     - sum(g_status[gr,j] for j=t:min.(t+generators_com[gr,:min_down_time]-1,T))
        #                     >= (- generators_com[gr,:min_down_time].*g_status[gr,t]
        #                     + generators_com[gr,:min_down_time].*generators_com[gr, :initial_status]))
        # [gr=G_com,t=2:T], (sum(g_status[gr,j] for j=t:min.(t+generators_com[gr,:min_down_time]-1,T))
        #                     >= - generators_com[gr,:min_down_time].*g_status[gr,t]
        #                     + generators_com[gr,:min_down_time].*g_status[gr,t-1])
    # end)



# 2. add all lines to the model
    # 2.1 set different lines types
    fix_lines_b = (.!network.lines[:s_nom_extendable])
    ext_lines_b = .!fix_lines_b
    lines = network.lines 

    # 2.2 iterator bounds
    LN_fix = sum(fix_lines_b)
    LN_ext = sum(ext_lines_b)

    # 2.3 add line variables to the model
    @variables m begin
        -lines[fix_lines_b,:s_nom][l]  <=  ln_fix[l=1:LN_fix,t=1:T] <= lines[fix_lines_b,:s_nom][l]
        ln_ext[l=1:LN_ext,t=1:T]
        lines[ext_lines_b,:s_nom_min][l] <=  ln_s_nom[l=1:LN_ext] <= lines[ext_lines_b,:s_nom_max][l]
    end

    # 2.4 add line constraint for extendable lines
    @constraints(m, begin
            [l=1:LN_ext,t=1:T], ln_ext[l,t] <=  ln_s_nom[l]
            [l=1:LN_ext,t=1:T], ln_ext[l,t] >= -ln_s_nom[l]
    end)

    ln = [ln_fix; ln_ext]
    lines = [lines[fix_lines_b,:]; lines[ext_lines_b,:]]


# 3. add all links to the model
    # 3.1 set different link types
    fix_links_b = .!network.links[:p_nom_extendable]
    ext_links_b = .!fix_links_b
    links = network.links

    # 3.2 iterator bounds
    LK_fix = sum(fix_links_b)
    LK_ext = sum(ext_links_b)

    #  3.3 set link variables
    @variables m begin
       ((links[fix_links_b, :p_nom].*links[fix_links_b, :p_min_pu])[l]  <=  lk_fix[l=1:LK_fix,t=1:T]
                <= (links[fix_links_b, :p_nom].*links[fix_links_b, :p_max_pu])[l])
        lk_ext[l=1:LK_ext,t=1:T]
        links[ext_links_b, :p_nom_min][l] <=  lk_p_nom[l=1:LK_ext] <= links[ext_links_b, :p_nom_max][l]
    end
    
    # 3.4 set constraints for extendable links
    @constraints(m, begin
    [l=1:LK_ext,t=1:T], lk_ext[l,t] >= lk_p_nom[l].*links[ext_links_b, :p_min_pu][l]
    [l=1:LK_ext,t=1:T], lk_ext[l,t] <= lk_p_nom[l].*links[ext_links_b, :p_max_pu][l]
    end)

    lk = [lk_fix; lk_ext]
    links = [links[fix_links_b,:]; links[ext_links_b,:]]

# 4. define storage_units
    # 4.1 set different storage_units types
    fix_su_b = .!network.storage_units[:p_nom_extendable]
    ext_su_b = .!fix_su_b
        # storage_units_spill = network.storage_units[inflow.max()>0,:]
    storage_units = network.storage_units

    # 4.2 iterator bounds
    SU_fix = sum(fix_su_b)
    SU_ext = sum(ext_su_b)
    # SU_spill = nrow(storage_units_spill)
    SU = nrow(storage_units)

    #  4.3 set link variables
    @variables m begin
       (0 <=  su_dispatch_fix[s=1:SU_fix,t=1:T] <=
                (storage_units[fix_su_b, :p_nom].*storage_units[fix_su_b, :p_max_pu])[s])
        su_dispatch_ext[s=1:SU_ext,t=1:T] >= 0
        (0 <=  su_store_fix[s=1:SU_fix,t=1:T] <=
                 - (storage_units[fix_su_b, :p_nom].*storage_units[fix_su_b, :p_min_pu])[s])
        su_store_ext[s=1:SU_ext,t=1:T] >= 0

        su_p_nom[s=1:SU_ext] >= 0

        0 <= su_soc_fix[s=1:SU_fix,t=1:T] <= (storage_units[fix_su_b,:max_hours]
                                            .*storage_units[fix_su_b,:p_nom])[s]
        su_soc_ext[s=1:SU_ext,t=1:T] >= 0

        # 0 <=  su_spill[l=1:SU_spill,t=1:T] <= inflow[l=1:SU_spill,t=1:T]
        end



    # 4.4 set constraints for extendable storage_units

    @constraints(m, begin
            [s=1:SU_ext,t=1:T], su_dispatch_ext[s,t] <= su_p_nom[s].*storage_units[ext_su_b, :p_max_pu][s]
            [s=1:SU_ext,t=1:T], su_store_ext[s,t] <= - su_p_nom[s].*storage_units[ext_su_b, :p_min_pu][s]
            [s=1:SU_ext,t=1:T], su_soc_ext[s,t] <= su_p_nom[s].*storage_units[ext_su_b, :max_hours][s]
    end)

    # 4.5 set charging constraint
    su_dispatch = [su_dispatch_fix; su_dispatch_ext]
    su_store = [su_store_fix; su_store_ext]
    su_soc = [su_soc_fix; su_soc_ext]
    storage_units = [storage_units[fix_su_b,:]; storage_units[ext_su_b,:]]  
    append_idx_col!(storage_units)

    is_cyclic_i = storage_units[storage_units[:cyclic_state_of_charge], :idx]
    not_cyclic_i = storage_units[.!storage_units[:cyclic_state_of_charge], :idx]

    @constraints(m, begin
            [s=is_cyclic_i,t=1], su_soc[s,t] == (su_soc[s,T]
                                        + storage_units[s,:efficiency_store] * su_store[s,t]
                                        - storage_units[s,:efficiency_dispatch] * su_dispatch[s,t])
            [s=not_cyclic_i,t=1], su_soc[s,t] == (storage_units[s,:state_of_charge_initial]
                                        + storage_units[s,:efficiency_store] * su_store[s,t]
                                        - storage_units[s,:efficiency_dispatch] * su_dispatch[s,t])

            [s=is_cyclic_i,t=2:T], su_soc[s,t] == (su_soc[s,t-1]
                                            + storage_units[s,:efficiency_store] * su_store[s,t]
                                            - storage_units[s,:efficiency_dispatch] * su_dispatch[s,t])

        end)


# # 5. define storages
#     # 5.1 set different storages types
#     fix_stores_b = .!network.storages[:e_nom_extendable]
#     ext_stores_b = .!fix_stores_b
#         # storages_spill = network.storages[inflow.max()>0,:]
#     storages = network.storages
#     append_idx_col!(storages)

#     # 5.2 iterator bounds
#     ST_fix = sum(fix_stores_b)
#     ST_ext = sum(ext_stores_b)
#     # ST_spill = nrow(storages_spill)
#     ST = nrow(storages)


#     #  5.3 set link variables
#     @variables m begin
#        (0 <=  st_dispatch_fix[s=1:ST_fix,t=1:T] <=
#                 (storages[fix_stores_b, :e_nom].*storages[fix_stores_b, :e_max_pu])[s])
#         st_dispatch_ext[s=1:ST_ext,t=1:T] >= 0
#         (0 <=  st_store_fix[s=1:ST_fix,t=1:T] <=
#                  - (storages[fix_stores_b, :e_nom].*storages[fix_stores_b, :e_min_pu])[s])
#         st_store_ext[s=1:ST_ext,t=1:T] >= 0

#         st_e_nom[s=1:ST_ext] >= 0

#         0 <= st_soc_fix[s=1:ST_fix,t=1:T] <= (storages[fix_stores_b,:max_hours]
#                                             .*storages[fix_stores_b,:e_nom])[s]
#         st_soc_ext[s=1:ST_ext,t=1:T] >= 0

#         # 0 <=  st_spill[l=1:ST_spill,t=1:T] <= inflow[l=1:ST_spill,t=1:T]
#         end



#     # 5.4 set constraints for extendable storages

#     @constraints(m, begin
#             [s=1:ST_ext,t=1:T], st_dispatch_ext[s,t] <= st_e_nom[s].*storages[ext_stores_b, :e_max_pu][s]
#             [s=1:ST_ext,t=1:T], st_store_ext[s,t] <= - st_e_nom[s].*storages[ext_stores_b, :e_min_pu][s]
#             [s=1:ST_ext,t=1:T], st_soc_ext[s,t] <= st_e_nom[s].*storages[ext_stores_b, :max_hours][s]
#     end)

#     # 5.5 set charging constraint
#     st_dispatch = [st_dispatch_fix; st_dispatch_ext]
#     st_store = [st_store_fix; st_store_ext]
#     st_soc = [st_soc_fix; st_soc_ext]
#     storages = [storages[fix_stores_b,:]; storages[ext_stores_b,:]]  

#     is_cyclic_i = storages[storages[:cyclic_state_of_charge], :idx]
#     not_cyclic_i = storages[.!storages[:cyclic_state_of_charge], :idx]

#     @constraints(m, begin
#             [s=is_cyclic_i,t=1], st_soc[s,t] == (st_soc[s,T]
#                                         + storages[s,:efficiency_store] * st_store[s,t]
#                                         - storages[s,:efficiency_dispatch] * st_dispatch[s,t])
#             [s=not_cyclic_i,t=1], st_soc[s,t] == (storages[s,:state_of_charge_initial]
#                                         + storages[s,:efficiency_store] * st_store[s,t]
#                                         - storages[s,:efficiency_dispatch] * st_dispatch[s,t])

#             [s=is_cyclic_i,t=2:T], st_soc[s,t] == (st_soc[s,t-1]
#                                             + storages[s,:efficiency_store] * st_store[s,t]
#                                             - storages[s,:efficiency_dispatch] * st_dispatch[s,t])

#         end)



## 6. define nodal balance constraint
    @constraint(m, balance[n=1:N, t=1:T], (
          sum(gn[findin(generators[:bus], [reverse_busidx[n]]), t])
        # + sum(gcom[idx_by(generators_com, :bus, [reverse_busidx[1]]), t])
        + sum(ln[ findin(lines[:bus1], [reverse_busidx[n]]) ,t])
        + sum(lk[ findin(links[:bus1], [reverse_busidx[n]]) ,t]) # *efficiency
        + sum(su_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

        - network.loads_t["p"][t,Symbol(reverse_busidx[n])]
        - sum(ln[ findin(lines[:bus0], [reverse_busidx[n]]) ,t])
        - sum(lk[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
        - sum(su_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

          == 0 ))

# 7. set Kirchhoff Voltage Law constraint
# since cyclebasis is not yet supported in LightGraphs, use the pyimported
# netwrokx in order to define all cycles. The cycle_basis returns a list of
# cycles, indicating the connected buses. For each cycle the connecting branches
# and their directions (signs) have to be determined.

# Might be nessecary to loop over all subgraphs as
# for (sn, sub) in enumerate(weakly_connected_components(g))
#     g_sub = induced_subgraph(g, sub)[1]

    append_idx_col!(lines)
    (branches, var, attribute) = (lines, ln, :x)
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
                        ln[cycles_branch[c],t]) == 0)
        # elseif attribute==:r
        #     @constraint(m, link_cycle_constraint[c=1:length(cycles_branch), t=1:T] ,
        #             dot(directions[c] .* links[cycles_branch[c], :r]/380. , lk[cycles_branch[c],t]) == 0)
        end
    end

# 8. set global_constraints
# only for co2_emissions till now

    if nrow(network.global_constraints)>0 && in("primary_energy", network.global_constraints[:_type])
        co2_limit = network.global_constraints[network.global_constraints[:name].=="co2_limit", :constant]
        nonnull_carriers = network.carriers[network.carriers[:co2_emissions].!=0, :]
        emmssions = Dict(zip(nonnull_carriers[:name], nonnull_carriers[:co2_emissions]))        
        carrier_index(carrier) = findin(generators[:carrier], [carrier])
        @constraint(m, sum(sum(dot(1./generators[carrier_index(carrier) , :efficiency], 
                    gn[carrier_index(carrier),t]) for t=1:T)  
                    * select_names(network.carriers, [carrier])[:co2_emissions]
                    for carrier in network.carriers[:name]) .<=  co2_limit)
    end

# 9. set objective function
    @objective(m, Min, 
                        sum(dot(generators[:marginal_cost], gn[:,t]) for t=1:T)
                        # consider already build infrastructur
                        + dot(generators[generators[:p_nom_extendable],:capital_cost], gen_p_nom[:]  )

                        + dot(lines[lines[:s_nom_extendable],:capital_cost], ln_s_nom[:])
                        + dot(links[links[:p_nom_extendable],:capital_cost], lk_p_nom[:])

                        + sum(dot(storage_units[:marginal_cost], su_dispatch[:,t]) for t=1:T)
                        + dot(storage_units[storage_units[:p_nom_extendable], :capital_cost], su_p_nom[:])
                        )

    status = solve(m)
# 10. extract optimisation results
    if status==:Optimal
        orig_gen_order = network.generators[:name]
        generators[generators[:p_nom_extendable],:p_nom] = getvalue(gen_p_nom)
        network.generators = generators
        network.generators_t["p"] = names!(DataFrame(transpose(getvalue(gn))),
                            Symbol.(generators[:name]))
        network.generators_t["p"] = network.generators_t["p"][:, Symbol.(orig_gen_order)]
        network.generators = select_names(network.generators, orig_gen_order)


        orig_line_order = network.lines[:name]
        network.lines = lines
        network.lines[BitArray(lines[:s_nom_extendable]),:s_nom] = getvalue(ln_s_nom)
        network.lines_t["p0"] = names!(DataFrame(transpose(getvalue(ln))),
                            Symbol.(lines[:name]))
        network.lines_t["p0"] = network.lines_t["p0"][:, Symbol.(orig_line_order)]
        network.lines = select_names(network.lines, orig_line_order)

        # network.buses_t["p"] =  DataFrame(ncols=nrow(network.buses))

        if nrow(links)>0
            orig_link_order = network.links[:name]
            network.links = links
            network.links[BitArray(links[:p_nom_extendable]),:p_nom] = getvalue(lk_p_nom)
            network.links_t["p0"] = names!(DataFrame(transpose(getvalue(lk))),
                                Symbol.(links[:name]))
            network.links_t["p0"] = network.links_t["p0"][:, Symbol.(orig_link_order)]
            network.links = select_names(network.links, orig_link_order)

        end
        if nrow(storage_units)>0
            orig_su_order = network.storage_units[:name]
            network.storage_units = storage_units
            network.storage_units[BitArray(storage_units[:p_nom_extendable]),:p_nom] = getvalue(su_p_nom)
            network.storage_units_t["p"] = names!(DataFrame(transpose(getvalue(su_dispatch .- su_store))),
                                Symbol.(storage_units[:name]))
            network.storage_units_t["p0"] = network.storage_units_t["p0"][:, Symbol.(orig_su_order)]
            network.storage_units = select_names(network.storage_units, orig_su_order)
        end
        println("Reduce cost to $(m.objVal)")
    end
    return m
end
