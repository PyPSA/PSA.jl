using JuMP
using MathProgBase

function run_lopf(network, solver; formulation::String="angles", objective::String="total", investment_type::String="continuous", blockmodel::Bool=false, decomposition::String="")

    if blockmodel
        println("Build block JuMP model.")
        m = build_block_lopf(network, solver; formulation=formulation, objective=objective, investment_type=investment_type,decomposition=decomposition)
    else
        println("Build ordinary JuMP model.")
        m = build_lopf(network, solver; formulation=formulation, objective=objective, investment_type=investment_type)
    end

    status = solve(m)

    if status==:Optimal

        # TODO: may be unneccesary (either here or in build_lopf)
        buses = network.buses
        lines = network.lines
        generators = network.generators
        links = network.links
        stores = network.stores
        storage_units = network.storage_units
        N = nrow(network.buses)
        L = nrow(network.lines)
        T = nrow(network.snapshots)
        fix_gens_b = (.!generators[:p_nom_extendable])
        ext_gens_b = convert(BitArray, generators[:p_nom_extendable])
        fix_lines_b = (.!lines[:s_nom_extendable])
        ext_lines_b = .!fix_lines_b
        fix_links_b = .!links[:p_nom_extendable]
        ext_links_b = .!fix_links_b
        fix_sus_b = .!storage_units[:p_nom_extendable]
        ext_sus_b = .!fix_sus_b
        fix_stores_b = .!stores[:e_nom_extendable]
        ext_stores_b = .!fix_stores_b

        G = [m[:G_fix]; m[:G_ext]]
        LN = [m[:LN_fix]; m[:LN_ext]]
        LK = [m[:LK_fix]; m[:LK_ext]]
        SU_dispatch = [m[:SU_dispatch_fix]; m[:SU_dispatch_ext]]
        SU_store = [m[:SU_store_fix]; m[:SU_store_ext]]
        SU_soc = [m[:SU_soc_fix]; m[:SU_soc_ext]]
        SU_spill = [m[:SU_spill_fix]; m[:SU_spill_ext]]
        ST_dispatch = [m[:ST_dispatch_fix]; m[:ST_dispatch_ext]]
        ST_store = [m[:ST_store_fix]; m[:ST_store_ext]]
        ST_soc = [m[:ST_soc_fix]; m[:ST_soc_ext]]
        ST_spill = [m[:ST_spill_fix], m[:ST_spill_ext]]

        lines = [lines[fix_lines_b,:]; lines[ext_lines_b,:]]
        generators = [generators[fix_gens_b,:]; generators[ext_gens_b,:] ]
        links = [links[fix_links_b,:]; links[ext_links_b,:]]
        storage_units = [storage_units[fix_sus_b,:]; storage_units[ext_sus_b,:]]
        stores = [stores[fix_stores_b,:]; stores[ext_stores_b,:]]

        println("The cost of transmission network are ", dot(lines[ext_lines_b,:capital_cost], getvalue(m[:LN_s_nom])) + dot(lines[fix_lines_b,:capital_cost], lines[fix_lines_b,:s_nom]))

        orig_gen_order = network.generators[:name]
        generators[:p_nom_opt] = deepcopy(generators[:p_nom])
        generators[ext_gens_b,:p_nom_opt] = getvalue(m[:G_p_nom])
        network.generators = generators
        network.generators_t["p"] = names!(DataFrame(transpose(getvalue(G))), Symbol.(generators[:name]))
        network.generators = select_names(network.generators, orig_gen_order)

        orig_line_order = network.lines[:name]
        network.lines = lines
        lines[:s_nom_opt] = deepcopy(lines[:s_nom])
        network.lines[ext_lines_b,:s_nom_opt] = getvalue(m[:LN_s_nom])
        network.lines_t["p0"] = names!(DataFrame(transpose(getvalue(LN))), Symbol.(lines[:name]))
        network.lines = select_names(network.lines, orig_line_order)

        # network.buses_t["p"] =  DataFrame(ncols=nrow(network.buses))

        if nrow(links)>0
            orig_link_order = network.links[:name]
            network.links = links
            links[:p_nom_opt] = deepcopy(links[:p_nom])
            network.links[ext_links_b,:p_nom_opt] = getvalue(m[:LK_p_nom])
            network.links_t["p0"] = names!(DataFrame(transpose(getvalue(LK))), Symbol.(links[:name]))
            network.links = select_names(network.links, orig_link_order)

        end
        if nrow(storage_units)>0
            orig_sus_order = network.storage_units[:name]
            network.storage_units = storage_units

            storage_units[:p_nom_opt] = deepcopy(storage_units[:p_nom])
            network.storage_units[ext_sus_b,:p_nom_opt] = getvalue(m[:SU_p_nom])
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