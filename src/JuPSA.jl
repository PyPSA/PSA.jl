module JuPSA

using DataFrames, JuMP, LightGraphs, Query, Gurobi, GraphLayout

using PyCall
const networkx = PyNULL()
copy!(networkx, pyimport("networkx" ))


export Network, lopf, import_network, idx, rev_idx, select_names, select_by, idx_by, to_symbol, append_idx_col!


mutable struct Network
    buses::DataFrame
    generators::DataFrame
    loads::DataFrame
    lines::DataFrame
    links::DataFrame
    storage_units::DataFrame
    transformers::DataFrame
    buses_t::Dict{String,DataFrame}
    generators_t::Dict{String,DataFrame}
    loads_t::Dict{String,DataFrame}
    lines_t::Dict{String,DataFrame}
    links_t::Dict{String,DataFrame}
    storage_units_t::Dict{String,DataFrame}
    transformers_t::Dict{String,DataFrame}
end



function Network(
    buses=DataFrame(repeat([Bool[]], outer=13),
[:name, :area, :area_offshore, :carrier, :control,
:country, :type, :v_nom, :v_mag_pu_set, :v_mag_pu_min,
:v_mag_pu_max, :x, :y]),

    generators=DataFrame(repeat([Bool[]], outer=27),
        [:bus, :control, :type, :p_nom, :p_nom_extendable,
        :p_nom_min, :p_nom_max, :p_min_pu, :p_max_pu, :p_set, :q_set, :sign,
        :carrier, :marginal_cost, :capital_cost, :efficiency, :committable,
        :start_up_cost, :shut_down_cost, :min_up_time, :min_down_time,
        :initial_status, :ramp_limit_up, :ramp_limit_down, :ramp_limit_start_up,
        :ramp_limit_shut_down, :p_nom_opt]),
    loads=DataFrame(repeat([Bool[]], outer=5),
        [:bus, :type, :p_set, :q_set, :sign]),
    lines=DataFrame(repeat([Bool[]], outer=23),
        [:bus0, :bus1, :type, :x, :r, :g, :b, :s_nom, :s_nom_extendable,
        :s_nom_min, :s_nom_max, :capital_cost, :length, :terrain_factor,
        :num_parallel, :v_ang_min, :v_ang_max, :sub_network, :x_pu,
        :r_pu, :g_pu, :b_pu, :s_nom_opt]) ,

    links=DataFrame(repeat([Bool[]], outer=16),
        [:bus0, :bus1, :type, :efficiency, :p_nom,
        :p_nom_extendable, :p_nom_min, :p_nom_max, :p_set, :p_min_pu,
        :p_max_pu, :capital_cost, :marginal_cost, :length, :terrain_factor,
        :p_nom_opt]),

    storage_units=DataFrame(repeat([Bool[]], outer=24),
        [:bus, :control, :type, :p_nom, :p_nom_extendable,
        :p_nom_min, :p_nom_max, :p_min_pu, :p_max_pu, :p_set,
        :q_set, :sign, :carrier, :marginal_cost, :capital_cost,
        :state_of_charge_initial, :state_of_charge_set,
        :cyclic_state_of_charge, :max_hours, :efficiency_store,
        :efficiency_dispatch, :standing_loss, :inflow, :p_nom_opt]),

    transformers=DataFrame(repeat([Bool[]], outer=26),
        [:bus0, :bus1, :type, :model, :x, :r, :g, :b, :s_nom,
        :s_nom_extendable, :s_nom_min, :s_nom_max, :capital_cost,
        :num_parallel, :tap_ratio, :tap_side, :tap_position,
        :phase_shift, :v_ang_min, :v_ang_max, :sub_network,
        :x_pu, :r_pu, :g_pu, :b_pu, :s_nom_opt]),

    buses_t=Dict([("marginal_price",DataFrame()), ("v_ang", DataFrame()),
            ("v_mag_pu_set",DataFrame()), ("q", DataFrame()),
            ("v_mag_pu", DataFrame()), ("p", DataFrame())]),
    generators_t=Dict{String,DataFrame}(
            [("marginal_price",DataFrame()), ("v_ang", DataFrame()),
            ("v_mag_pu_set",DataFrame()), ("q", DataFrame()),
            ("v_mag_pu", DataFrame()), ("p", DataFrame())]),
    loads_t=Dict{String,DataFrame}([("q_set",DataFrame()), ("p_set", DataFrame()),
            ("q", DataFrame()), ("p", DataFrame())]),
    lines_t=Dict{String,DataFrame}([("q0",DataFrame()), ("q1", DataFrame()),
            ("p0",DataFrame()), ("p1", DataFrame()),
            ("mu_lower", DataFrame()), ("mu_upper", DataFrame())]),
    links_t=Dict{String,DataFrame}([("p_min_pu",DataFrame()), ("p_max_pu", DataFrame()),
            ("p0",DataFrame()), ("p1", DataFrame()),
            ("mu_lower", DataFrame()), ("mu_upper", DataFrame()),
            ("efficiency",DataFrame()), ("p_set", DataFrame())]),
    storage_units_t=Dict{String,DataFrame}([("p_min_pu",DataFrame()), ("p_max_pu", DataFrame()),
        ("inflow",DataFrame()), ("mu_lower", DataFrame()), ("mu_upper", DataFrame()),
        ("efficiency",DataFrame()), ("p_set", DataFrame())]),
    transformers_t=Dict{String,DataFrame}([("p0",DataFrame()), ("p1", DataFrame()),
            ("q0", DataFrame()), ("q1", DataFrame()),
            ("mu_upper",DataFrame())])
    )
    Network(
        buses, generators, loads, lines, links, transformers, buses_t,
        generators_t, loads_t, lines_t, links_t, storage_unists_t, transformers_t)
end

function import_network(file="")
    n = Network()
    if file!=""
        if ispath("$file/buses.csv")
            n.buses = readtable("$file/buses.csv", truestrings=["True"],
                                            falsestrings=["False"])
        end
        if ispath("$file/lines.csv")
            n.lines = readtable("$file/lines.csv", truestrings=["True"],
                                            falsestrings=["False"])
        end
        if ispath("$file/links.csv")
            n.links = readtable("$file/links.csv", truestrings=["True"],
                                            falsestrings=["False"])
        end
        if ispath("$file/generators.csv")
            n.generators = readtable("$file/generators.csv", truestrings=["True"],
                                            falsestrings=["False"])
        end
        if ispath("$file/loads.csv")
            n.loads = readtable("$file/loads.csv", truestrings=["True"],
                                            falsestrings=["False"])
        end
        if ispath("$file/storage_units.csv")
            n.storage_units = readtable("$file/storage_units.csv", truestrings=["True"],
                                            falsestrings=["False"])
        end
        if ispath("$file/loads-p_set.csv")
            n.loads_t["p"] = readtable("$file/loads-p_set.csv")
        end
        if ispath("$file/generators-p_max_pu.csv")
            n.generators_t["p_max_pu"] = readtable("$file/generators-p_max_pu.csv")
        end
        if ispath("$file/storage_units-inflow.csv")
            n.storage_units_t["inflow"] = readtable("$file/storage_units-inflow.csv")
        end
    end
    return n
end



function set_snapshots(network, snapshots)
    for field=[field for field=fieldnames(network) if String(field)[end-1:end]=="_t"]
        for df_name=keys(getfield(network,field))
            if nrow(getfield(network,field)[df_name])>0
                getfield(network,field)[df_name] = getfield(network,field)[df_name][snapshots, :]
            end
        end
    end
end

# auxilliary functions
idx(dataframe) = Dict(zip(dataframe[:name], Iterators.countfrom(1)))
rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe[:name]))
select_names(dataframe, names) = dataframe[findin(dataframe[:name], names),:]
select_by(dataframe, col, values) = dataframe[findin(dataframe[col], values),:]
idx_by(dataframe, col, values) = select_by(dataframe, col, values)[:idx]
function to_symbol(str)
    if typeof(str)==String
        return Symbol(replace(str, " ", "_"))
    elseif typeof(str)==Vector{String}
        return [Symbol(replace(x, " ", "_")) for x in str]
    elseif typeof(str)==Int
        return Symbol("$str")
    else
        return str
    end
end

function to_string(sym)
    if typeof(sym)==Symbol
        return replace(String(sym), "_", " ")
    # else
    #     return sym
    end
end

function append_idx_col!(dataframe)
    if typeof(dataframe)==Vector{DataFrames.DataFrame}
        for df in dataframe
            df[:idx] = 1:nrow(df)
        end
    else
        dataframe[:idx] = 1:nrow(dataframe)
    end
end

function get_switchable_as_dense(network, component, attribute, snapshots=0)
    if snapshots==0
        snapshots = network.loads_t["p"][:name] # network.snapshots
    end
    T = length(snapshots)
    component_t = to_symbol(component * "_t")
    component = to_symbol(component)
    dense = DataFrame()
    if in(attribute, keys(getfield(network, component_t)))
        dense = getfield(network, component_t)[attribute]
    end
    cols = to_symbol(collect(getfield(network, component)[:name]))
    not_included = [to_string(c) for c=cols if !in(c,names(dense))]
    if length(not_included)>0
        attribute = to_symbol(attribute)
        df = select_names(getfield(network, component), not_included)
        df = names!(DataFrame(repmat(transpose(Array(df[attribute])), T)),
                to_symbol(not_included))
        dense = [dense df]
    end
    return dense[cols]
end




function calculate_dependent_values(network)
    function set_default(dataframe, col, default)
        if !in(col, names(dataframe))
            dataframe[col] = default
        end
    end
    defaults = [(:p_nom_extendable, false), (:p_nom_max, Inf),(:commitable, false),
                (:p_min_pu, 0), (:p_max_pu, 1), (:p_nom_min, 0),(:capital_cost, 0),
                (:min_up_time, 0), (:min_down_time, 0), (:initial_status, true)]
    for (col, default) in defaults
        set_default(network.generators, col, default)
    end
    defaults = [(:s_nom_extendable, false), (:s_nom_min, 0),(:s_nom_max, Inf),
                (:s_nom_min, 0), (:s_nom_max, Inf), (:capital_cost, 0)]
    for (col, default) in defaults
        set_default(network.lines, col, default)
    end
    defaults = [(:p_nom_extendable, false), (:p_nom_max, Inf), (:p_min_pu, 0),
                (:p_max_pu, 1),(:p_nom_min, 0), (:p_nom_max, Inf), (:capital_cost, 0)]
    for (col, default) in defaults
        set_default(network.links, col, default)
    end
        defaults = [(:p_nom_min, 0), (:p_nom_max, Inf), (:p_min_pu, -1),
                    (:p_max_pu, 1), (:marginal_cost, 0)]
        for (col, default) in defaults
            set_default(network.storage_units, col, default)
    end
end



function lopf(network)

    m = Model(solver=GurobiSolver())

    buses = network.buses
    reverse_busidx = rev_idx(buses)
    busidx = idx(buses)
    N = nrow(network.buses)
    T = nrow(network.loads_t["p"]) #normally snapshots


# 1. add all generators to the model
    # 1.1 set different generator types
    generators_fix = network.generators[(!network.generators[:p_nom_extendable]) .&
                                    (!network.generators[:commitable]),:]
    generators_ext = network.generators[network.generators[:p_nom_extendable],:]
    generators_com = network.generators[network.generators[:commitable],:]
    generators = [generators_fix; generators_ext; generators_com]
    append_idx_col!([generators_fix, generators_ext, generators_com, generators])

    # 1.2 fix bounds for iterating
    G_fix = nrow(generators_fix)
    G_ext = nrow(generators_ext)
    G_com = nrow(generators_com)

   # 1.3 add generator variables to the model
    @variables m begin

        (generators_fix[gr, :p_nom].*generators_fix[gr, :p_min_pu] <= g_fix[gr=1:G_fix,t=1:T]
                                <= generators_fix[gr, :p_nom].*generators_fix[gr, :p_max_pu])

        g_ext[gr=1:G_ext,t=1:T]
        generators_ext[gr, :p_nom_min] <=  gen_p_nom[gr=1:G_ext] <= generators_ext[gr, :p_nom_max]

        g_status[gr=1:G_com,t=1:T], Bin
        (generators_com[gr, :p_nom].*generators_com[gr, :p_min_pu] <= g_com[gr=1:G_com,t=1:T]
                                <= generators_com[gr, :p_nom].*generators_com[gr, :p_max_pu])
    end

    gn = [g_fix; g_ext; g_com] # gn is the concatenated variable array

    # 1.4 set constraints for generators

    g_up_time_i = generators_com[generators_com[:min_up_time].>0, :idx]
    g_down_time_i = generators_com[generators_com[:min_down_time].>0, :idx]

    @constraints(m, begin
        [gr=1:G_ext,t=1:T], g_ext[gr,t] >= gen_p_nom[gr].*generators_ext[gr, :p_min_pu]
        [gr=1:G_ext,t=1:T], g_ext[gr,t] <= gen_p_nom[gr].*generators_ext[gr, :p_max_pu]

        [gr=1:G_com,t=1:T], g_com[gr,t] - g_com[gr,t].*g_status[gr,t] == 0

        # [gr=g_up_time_i,t=1], (sum(g_status[gr,j] for j=t:min.(t+generators_com[gr,:min_up_time]-1,T))
        #                     >=
        #                     # generators_com[gr,:min_up_time].*g_status[gr,t]
        #                     generators_com[gr,:min_up_time].*generators_com[gr, :initial_status])
        # [gr=g_up_time_i,t=2:T], (sum(g_status[gr,j] for j=t:min.(t+generators_com[gr,:min_up_time]-1,T))
        #                     >= generators_com[gr,:min_up_time].*g_status[gr,t]
        #                     - generators_com[gr,:min_up_time].*g_status[gr,t-1])
        #
        # [gr=g_down_time_i,t=1], (generators_com[gr,:min_down_time]
        #                     - sum(g_status[gr,j] for j=t:min.(t+generators_com[gr,:min_down_time]-1,T))
        #                     >= (- generators_com[gr,:min_down_time].*g_status[gr,t]
        #                     + generators_com[gr,:min_down_time].*generators_com[gr, :initial_status]))
        # [gr=g_down_time_i,t=2:T], (sum(g_status[gr,j] for j=t:min.(t+generators_com[gr,:min_down_time]-1,T))
        #                     >= - generators_com[gr,:min_down_time].*g_status[gr,t]
        #                     + generators_com[gr,:min_down_time].*g_status[gr,t-1])
    end)


# 2. add all lines to the model
    # 2.1 set different lines types
    lines_fix = network.lines[!network.lines[:s_nom_extendable],:]
    lines_ext = network.lines[network.lines[:s_nom_extendable],:]
    lines = [lines_fix; lines_ext]
    append_idx_col!([lines_fix, lines_ext, lines])

    # 2.2 iterator bounds
    LN_fix = nrow(lines_fix)
    LN_ext = nrow(lines_ext)

    # 2.3 add line variables to the model
    @variables m begin
        -lines_fix[l,:s_nom]  <=  ln_fix[l=1:LN_fix,t=1:T] <= lines_fix[l,:s_nom]
        ln_ext[l=1:LN_ext,t=1:T]
        lines_ext[l,:s_nom_min] <=  ln_s_nom[l=1:LN_ext] <= lines_ext[l,:s_nom_max]
    end

    ln = [ln_fix; ln_ext]

    # 2.4 add line constraint for extendable lines
    @constraints(m, begin
            [l=1:LN_ext,t=1:T], ln_ext[l,t] <=  ln_s_nom[l]
            [l=1:LN_ext,t=1:T], ln_ext[l,t] >= -ln_s_nom[l]
    end)


# 3. add all links to the model
    # 3.1 set different link types
    links_fix = network.links[!network.links[:p_nom_extendable],:]
    links_ext = network.links[network.links[:p_nom_extendable],:]
    links = [links_fix; links_ext]
    append_idx_col!([links_fix, links_ext, links])

    # 3.2 iterator bounds
    LK_fix = nrow(links_fix)
    LK_ext = nrow(links_ext)

    #  3.3 set link variables
    @variables m begin
       (links_fix[l, :p_nom].*links_fix[l, :p_min_pu]  <=  lk_fix[l=1:LK_fix,t=1:T]
                <= links_fix[l, :p_nom].*links_fix[l, :p_max_pu])
        lk_ext[l=1:LK_ext,t=1:T]
        links_ext[l, :p_nom_min] <=  lk_p_nom[l=1:LK_ext] <= links_ext[l, :p_nom_max]
    end
    lk = [lk_fix; lk_ext]

    # 3.4 set constraints for extendable links
    @constraints(m, begin
            [l=1:LK_ext,t=1:T], lk_ext[l,t] >= lk_p_nom[l].*links_ext[l, :p_min_pu]
            [l=1:LK_ext,t=1:T], lk_ext[l,t] <= lk_p_nom[l].*links_ext[l, :p_max_pu]
    end)


# 4. define storage_units
    # 4.1 set different storage_units types
    storage_units_fix = network.storage_units[.!network.storage_units[:p_nom_extendable],:]
    storage_units_ext = network.storage_units[network.storage_units[:p_nom_extendable],:]
        # storage_units_spill = network.storage_units[inflow.max()>0,:]
    storage_units = [storage_units_fix; storage_units_ext]
    append_idx_col!([storage_units_fix, storage_units_ext, # storage_units_spill,
                    storage_units])

    # 4.2 iterator bounds
    SU_fix = nrow(storage_units_fix)
    SU_ext = nrow(storage_units_ext)
    # SU_spill = nrow(storage_units_spill)
    SU = nrow(storage_units)

    #  4.3 set link variables
    @variables m begin
       (0 <=  su_dispatch_fix[s=1:SU_fix,t=1:T] <=
                storage_units_fix[s, :p_nom].*storage_units_fix[s, :p_max_pu])
        su_dispatch_ext[s=1:SU_ext,t=1:T] >= 0
        (0 <=  su_store_fix[s=1:SU_fix,t=1:T] <=
                 - storage_units_fix[s, :p_nom].*storage_units_fix[s, :p_min_pu])
        su_store_ext[s=1:SU_ext,t=1:T] >= 0

        su_p_nom[s=1:SU_ext] >= 0

        (0 <= su_soc_fix[s=1:SU_fix,t=1:T] <= storage_units_fix[s,:max_hours]
                                            .*storage_units_fix[s,:p_nom])
        su_soc_ext[s=1:SU_ext,t=1:T] >= 0

        # 0 <=  su_spill[l=1:SU_spill,t=1:T] <= inflow[l=1:SU_spill,t=1:T]
        end
    su_dispatch = [su_dispatch_fix; su_dispatch_ext]
    su_store = [su_store_fix; su_store_ext]
    su_soc = [su_soc_fix; su_soc_ext]


    # 4.4 set constraints for extendable storage_units
    is_cyclic_i = storage_units[storage_units[:cyclic_state_of_charge], :idx]
    not_cyclic_i = storage_units[.!storage_units[:cyclic_state_of_charge], :idx]

    @constraints(m, begin
            [s=1:SU_ext,t=1:T], su_dispatch_ext[s,t] <= su_p_nom[s].*storage_units_ext[s, :p_max_pu]
            [s=1:SU_ext,t=1:T], su_store_ext[s,t] <= - su_p_nom[s].*storage_units_ext[s, :p_min_pu]
            [s=1:SU_ext,t=1:T], su_soc_ext[s,t] <= su_p_nom[s].*storage_units_ext[s, :max_hours]

            [s=is_cyclic_i,t=1], su_soc[s,t] == su_soc[s,T] + su_store[s,t] - su_dispatch[s,t]
            [s=not_cyclic_i,t=1], su_soc[s,t] == storage_units[s,:state_of_charge_initial] + su_store[s,t] - su_dispatch[s,t]

            [s=is_cyclic_i,t=2:T], su_soc[s,t] == su_soc[s,t-1] + su_store[s,t] - su_dispatch[s,t]
        end)


## 5. define nodal balance constraint
    @constraint(m, balance[n=1:N, t=1:T], (
          sum(g_fix[idx_by(generators_fix, :bus, [reverse_busidx[n]]), t])
        + sum(g_ext[idx_by(generators_ext, :bus, [reverse_busidx[n]]), t])
        # + sum(gcom[idx_by(generators_com, :bus, [reverse_busidx[1]]), t])
        + sum(ln[ idx_by(lines, :bus1, [reverse_busidx[n]]) ,t])
        + sum(lk[ idx_by(links, :bus1, [reverse_busidx[n]]) ,t]) # *efficiency
        + sum(su_dispatch[ idx_by(storage_units, :bus, [reverse_busidx[n]]) ,t])

        - network.loads_t["p"][t,to_symbol(reverse_busidx[n])]
        - sum(ln[ idx_by(lines, :bus0, [reverse_busidx[n]]) ,t])
        - sum(lk[ idx_by(links, :bus0, [reverse_busidx[n]]) ,t])
        - sum(su_store[ idx_by(storage_units, :bus, [reverse_busidx[n]]) ,t])

          == 0 ))

# 6. set Kirchhoff Voltage Law constraint
# since cyclebasis is not yet supported in LightGraphs, use the pyimported
# netwrokx in order to define all cycles. The cycle_basis returns a list of
# cycles, indicating the connected buses. For each cycle the connecting branches
# and their directions (signs) have to be determined.

# Might be nessecary to loop over all subgraphs as
# for (sn, sub) in enumerate(weakly_connected_components(g))
#     # g_sub = induced_subgraph(g, sub)[1]


    (branches, var, attribute) = (lines, ln, :x)
    g = networkx[:Graph]()
    g[:add_nodes_from](busidx)
    g[:add_edges_from]([(busidx[l[:bus0]], busidx[l[:bus1]]) for l in eachrow(branches)])
    cycles = networkx[:cycle_basis](g)
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
                    dot(directions[c] .* lines[cycles_branch[c], :x]/380., ln[cycles_branch[c],t]) == 0)
        # elseif attribute==:r
        #     @constraint(m, link_cycle_constraint[c=1:length(cycles_branch), t=1:T] ,
        #             dot(directions[c] .* links[cycles_branch[c], :r]/380. , lk[cycles_branch[c],t]) == 0)
        end
    end


# 7. set objective function
    @objective(m, Min, sum(dot(generators[:marginal_cost], gn[:,t]) for t=1:T)
                        + dot(generators_ext[:capital_cost], gen_p_nom[:])

                        + dot(lines_ext[:capital_cost], ln_s_nom[:])
                        + dot(links_ext[:capital_cost], lk_p_nom[:])

                        + sum(dot(storage_units[:marginal_cost], su_dispatch[:,t]) for t=1:T)
                        + dot(storage_units_ext[:capital_cost], su_p_nom[:])
                        )

    status = solve(m)
# 8. extract optimisation results
    if status==:Optimal
        network.generators[:p_nom] = DataArray{Float64}(network.generators[:p_nom])
        network.generators[network.generators[:p_nom_extendable],:p_nom] = getvalue(gen_p_nom)
        network.generators_t["p"] = names!(DataFrame(transpose(getvalue(gn))),
                            [to_symbol(n) for n=generators[:name]])

        network.lines[:s_nom] = DataArray{Float64}(network.lines[:s_nom])
        network.lines[network.lines[:s_nom_extendable],:s_nom] = getvalue(ln_s_nom)
        network.lines_t["p0"] = names!(DataFrame(transpose(getvalue(ln))),
                            [to_symbol(n) for n=lines[:name]])
        if length(eachrow(links))>0
            network.links[:p_nom] = DataArray{Float64}(network.links[:p_nom])
            network.links[network.links[:p_nom_extendable],:p_nom] = getvalue(lk_p_nom)
            network.links_t["p0"] = names!(DataFrame(transpose(getvalue(lk))),
                                [to_symbol(n) for n=links[:name]])
        end
        if length(eachrow(storage_units))>0
            network.storage_units[:p_nom] = DataArray{Float64}(network.storage_units[:p_nom])
            network.storage_units[network.storage_units[:p_nom_extendable],:p_nom] = getvalue(su_p_nom)
            network.storage_units_t["p"] = names!(DataFrame(transpose(getvalue(su_dispatch .- su_store))),
                                [to_symbol(n) for n=storage_units[:name]])
        end
    end

end

end
