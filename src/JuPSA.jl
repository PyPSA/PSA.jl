module JuPSA

using DataFrames, JuMP, LightGraphs, Query, Gurobi, GraphLayout
export Network, lopf, import_network

mutable struct Network
    buses::DataFrame
    generators::DataFrame
    loads::DataFrame
    lines::DataFrame
    links::DataFrame
    transformers::DataFrame
    buses_t::Dict{String,DataFrame}
    generators_t::Dict{String,DataFrame}
    loads_t::Dict{String,DataFrame}
    lines_t::Dict{String,DataFrame}
    links_t::Dict{String,DataFrame}
    transformers_t::Dict{String,DataFrame}
end

function Network(
    buses=DataFrame(repeat([Int64[]], outer=10), [:name, :area, :area_offshore,
        :carrier, :control, :country, :_type, :v_nom, :x, :y]),
    generators=DataFrame(repeat([Int64[]], outer=12), [:name, :bus, :capital_cost,
        :carrier, :control, :efficiency, :marginal_cost, :p_nom, :p_nom_extendable,
        :p_nom_max, :_type, :weight]),
    loads=DataFrame(repeat([Int64[]], outer=3), [:name, :bus, :_type]),
    lines=DataFrame(repeat([Int64[]], outer=10), [:name, :bus0, :bus1, :b, :length,
         :r, :s_nom, :_type, :x, :extendable]) ,
    links=DataFrame(repeat([Int64[]], outer=12), [:name, :bus0, :bus1, :circuits,
        :geometry, :length, :p_nom, :tags, :_type, :under_construction,
        :underground, :extendable]),
    transformers=DataFrame(),
    buses_t=Dict{String,DataFrame}(),
    generators_t=Dict{String,DataFrame}(),
    loads_t=Dict{String,DataFrame}(),
    lines_t=Dict{String,DataFrame}(),
    links_t=Dict{String,DataFrame}(),
    transformers_t=Dict{String,DataFrame}()
    )
    Network(
        buses,
        generators,
        loads,
        lines,
        links,
        transformers,
        buses_t,
        generators_t,
        loads_t,
        lines,
        links_t,
        transformers_t
    )
end

function import_network(file="")
    n = Network()
    if file!=""
        if ispath("$file/buses.csv")
            n.buses = readtable("$file/buses.csv")
        end
        if ispath("$file/lines.csv")
            n.lines = readtable("$file/lines.csv")
        end
        if ispath("$file/links.csv")
            n.links = readtable("$file/links.csv")
        end
        if ispath("$file/generators.csv")
            n.generators = readtable("$file/generators.csv")
        else
            n.generators = DataFrame(repeat([[NaN]], outer=9), [:name, :bus0, :bus1, :b, :length,
                 :r, :s_nom, :_type, :x])
        end
        if ispath("$file/loads.csv")
            n.loads = readtable("$file/loads.csv")
        end
        if ispath("$file/loads-p_set.csv")
            n.loads_t["p"] = readtable("$file/loads-p_set.csv")
        end
    end
    return n
end


# function set_snapshots(n::Network, Index)
#     for timedependent in [n.buses_t, n.lines_t, n.generators_t]
#         timedependent = loc(timedependent)[Index]
#     end
# end

idx(dataframe) = Dict(zip(dataframe[:name], Iterators.countfrom(1)))
rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe[:name]))

index_name(dataframe, names) = dataframe[findin(dataframe[:name], names),:]
index_by(dataframe, col, values) = dataframe[findin(dataframe[col], values),:]
idx_by(dataframe, col, values) = index_by(dataframe, col, values)[:idx]




function calculate_dependent_values(network)
    function set_default(dataframe, col, default)
        if !in(col, names(dataframe))
            dataframe[col] = default
        end
    end
    defaults = [(:extendable, false), (:commitable, false), (:p_min_pu, 0), (:p_max_pu, 1)]
    for (col, default) in defaults
        set_default(network.generators, col, default)
    end
    defaults = [(:extendable, false), (:s_nom_min, NaN),(:s_nom_max, NaN)]
    for (col, default) in defaults
        set_default(network.lines, col, default)
    end
    defaults = [(:extendable, false), (:p_min_pu, 0), (:p_max_pu, 1)]
    for (col, default) in defaults
        set_default(network.links, col, default)
    end
end



function lopf(network)
    buses = network.buses
    generators = network.generators[(!network.generators[:extendable]) .&
                                    (!network.generators[:commitable]),:]
    generators_ext = network.generators[network.generators[:extendable],:]
    generators_com = network.generators[network.generators[:commitable],:]
    lines = network.lines[!network.lines[:extendable],:]
    lines_ext = network.lines[network.lines[:extendable],:]
    links = network.links[!network.links[:extendable],:]
    links_ext = network.links[network.links[:extendable],:]


    for df in [buses, generators, generators_ext, generators_com, lines, lines_ext,
        links, links_ext]
        df[:idx] = 1:nrow(df)
    end

    N = nrow(network.buses)
    T = nrow(network.loads_t["p"])
    LN = nrow(lines)
    LN_ext = nrow(lines_ext)
    G = nrow(generators)
    G_ext = nrow(generators_ext)
    G_com = nrow(generators_com)
    LK = nrow(links)
    LK_ext = nrow(links_ext)


    reverse_busidx = rev_idx(buses)
    busidx = idx(buses)

    m = Model(solver=GurobiSolver())

    @variables m begin

        generators[gr, :p_nom].*generators[gr, :p_min_pu] <= gt[gr=1:G,t=1:T] <= generators[gr, :p_nom].*generators[gr, :p_max_pu]
        gext[gr=1:G_ext,t=1:T]
        generators_ext[gr, :p_nom_min]          <=  g_p_nom[gr=1:G_ext]      <= generators_ext[gr, :p_nom_max]

        -lines[l, :s_nom]                       <=  lnt[l=1:LN,t=1:T]        <= lines[l, :s_nom]
        lnext[l=1:LN_ext,t=1:T]
        lines_ext[l, :s_nom_min]                <=  ln_s_nom[l=1:LN_ext]     <= lines_ext[l, :s_nom_max]


        -links[l, :p_nom].*links[l, :p_min_pu]  <=  lkt[l=1:LK,t=1:T]        <= links[l, :p_nom].*links[l, :p_max_pu]
        lkext[l=1:LK_ext,t=1:T]
        links_ext[l, :p_nom_min]                <=  lk_p_nom[l=1:LK_ext]     <= links_ext[l, :p_nom_max]

    end

    @constraints(m, begin
            [gr=1:G_ext,t=1:T], gext[gr,t] >= g_p_nom[gr].*generators_ext[gr, :p_min_pu]
            [gr=1:G_ext,t=1:T], gext[gr,t] <= g_p_nom[gr].*generators_ext[gr, :p_max_pu]

            [l=1:LN_ext,t=1:T], lnext[l,t] <=  ln_s_nom[l]
            [l=1:LN_ext,t=1:T], lnext[l,t] >= -ln_s_nom[l]

            [l=1:LK_ext,t=1:T], lkext[l,t] >= lk_p_nom[l].*links_ext[l, :p_min_pu]
            [l=1:LK_ext,t=1:T], lkext[l,t] <= lk_p_nom[l].*links_ext[l, :p_max_pu]
    end)


    to_symbol(str) = Symbol(replace(str, " ", "_"))
    #   nodal balance
    @constraint(m, balance[n=1:N, t=1:T], (
          sum(gt[idx_by(generators, :bus, [reverse_busidx[n]]), t])
        + sum(gext[idx_by(generators_ext, :bus, [reverse_busidx[n]]), t])
        # + sum(gcom[idx_by(generators_com, :bus, [reverse_busidx[1]]), t])

        - network.loads_t["p"][t,to_symbol(reverse_busidx[n])]

        + sum(lnt[ idx_by(lines, :bus1, [reverse_busidx[n]]) ,t])
        - sum(lnt[ idx_by(lines, :bus0, [reverse_busidx[n]]) ,t])
        + sum(lnext[ idx_by(lines_ext, :bus1, [reverse_busidx[n]]) ,t])
        - sum(lnext[ idx_by(lines_ext, :bus0, [reverse_busidx[n]]) ,t])

        + sum(lkt[ idx_by(links, :bus1, [reverse_busidx[n]]) ,t])
        - sum(lkt[ idx_by(links, :bus0, [reverse_busidx[n]]) ,t])
        + sum(lkext[ idx_by(links_ext, :bus1, [reverse_busidx[n]]) ,t])
        - sum(lkext[ idx_by(links_ext, :bus0, [reverse_busidx[n]]) ,t])
          == 0 ))

    lnall=[lnt; lnext]
    linesall = [lines; lines_ext]

    g = DiGraph(length(busidx))
    for l = eachrow(linesall)
      add_edge!(g, busidx[l[:bus0]], busidx[l[:bus1]])
    end

    #   Kirchhoff voltage law
    cycles = [i for i in simplecycles(g) if length(i)>2]
    @constraint(m, cycle_constraint[c=1:length(cycles), t=1:T] ,
                sum(dot(linesall[cycles[c], :x]/380. , lnall[cycles[c],t])) == 0)

    # Might be nessecary to loop over all subgraphs as
    # for (sn, sub) in enumerate(weakly_connected_components(g))
        # g_sub = induced_subgraph(g, sub)[1]



    @objective(m, Min, sum(dot(generators[:marginal_cost], gt[:,t]) for t=1:T))

    slove(m)
    return m
end

end
