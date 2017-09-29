using DataFrames, JuMP, LightGraphs, Query, Gurobi, GraphLayout

using JuPSA





network = import_network("mininetwork")
network.loads_t["p"] = DataFrame(My_bus_0 = 0, My_bus_1 = 0, My_bus_2 = 100)

# network = JuPSA.import_network("pre8-37")
# network.loads_t["p"][:LT0_0] =0
# network.generators = network.generators[network.generators[:carrier].=="OCGT",:]
# network.generators[:p_nom] = network.generators[:p_nom] *100
# network.loads_t["p"] = network.loads_t["p"][1:1000,:]


idx(dataframe) = Dict(zip(dataframe[:name], Iterators.countfrom(1)))
rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe[:name]))

index_name(dataframe, names) = dataframe[findin(dataframe[:name], names),:]
index_by(dataframe, col, values) = dataframe[findin(dataframe[col], values),:]
idx_by(dataframe, col, values) = index_by(dataframe, col, values)[:idx]


JuPSA.calculate_dependent_values(network)
# network.lines[1, :extendable] = true
# network.lines[1, :s_nom_min] = 0
# network.lines[1, :s_nom_max] = 1000


buses = network.buses

if length(network.generators)>0
    generators_fix = network.generators[(!network.generators[:extendable]) .&
                                    (!network.generators[:commitable]),:]
    generators_ext = network.generators[network.generators[:extendable],:]
    generators_com = network.generators[network.generators[:commitable],:]
    generators = [generators_fix; generators_ext; generators_com]
else
    for df in [generators_fix, generators_ext, generators_com, generators]
        df = network.generators
    end
end

if length(network.lines)>0
    lines_fix = network.lines[!network.lines[:extendable],:]
    lines_ext = network.lines[network.lines[:extendable],:]
    lines = [lines_fix; lines_ext]
else
    lines_fix = network.lines
    lines_ext = network.lines
    lines = [lines_fix; lines_ext]
end

if length(network.links)>0
    links_fix = network.links[network.links[:extendable],:]
    links_ext = network.links[network.links[:extendable],:]
    links = [links_fix; links_ext]
else
    links_fix = network.links
    links_ext = network.links
    links = [links_fix; links_ext]
end

for df in [buses, generators_fix, generators_ext, generators_com, lines_fix, lines_ext,
    links_fix, links_ext]
    df[:idx] = 1:nrow(df)
end

N = nrow(network.buses)
T = nrow(network.loads_t["p"])
G_fix = nrow(generators_fix)
G_ext = nrow(generators_ext)
G_com = nrow(generators_com)
LN_fix = nrow(lines_fix)
LN_ext = nrow(lines_ext)
LK_fix = nrow(links_fix)
LK_ext = nrow(links_ext)



reverse_busidx = rev_idx(buses)
busidx = idx(buses)

m = Model(solver=GurobiSolver())

@variables m begin

    (generators_fix[gr, :p_nom].*generators_fix[gr, :p_min_pu] <= g_fix[gr=1:G_fix,t=1:T]
                            <= generators_fix[gr, :p_nom].*generators_fix[gr, :p_max_pu])
    g_ext[gr=1:G_ext,t=1:T]
    generators_ext[gr, :p_nom_min]          <=  g_p_nom[gr=1:G_ext]      <= generators_ext[gr, :p_nom_max]

    -lines_fix[l, :s_nom]                       <=  ln_fix[l=1:LN_fix,t=1:T]        <= lines_fix[l, :s_nom]
    ln_ext[l=1:LN_ext,t=1:T]
    lines_ext[l, :s_nom_min]                <=  ln_s_nom[l=1:LN_ext]     <= lines_ext[l, :s_nom_max]


    -links_fix[l, :p_nom].*links_fix[l, :p_min_pu]  <=  lk_fix[l=1:LK_fix,t=1:T]        <= links_fix[l, :p_nom].*links_fix[l, :p_max_pu]
    lk_ext[l=1:LK_ext,t=1:T]
    links_ext[l, :p_nom_min]                <=  lk_p_nom[l=1:LK_ext]     <= links_ext[l, :p_nom_max]

end

ln = [ln_fix; ln_ext]
lk = [lk_fix; lk_ext]
gn = [g_fix; g_ext]

@constraints(m, begin
        [gr=1:G_ext,t=1:T], g_ext[gr,t] >= g_p_nom[gr].*generators_ext[gr, :p_min_pu]
        [gr=1:G_ext,t=1:T], g_ext[gr,t] <= g_p_nom[gr].*generators_ext[gr, :p_max_pu]

        [l=1:LN_ext,t=1:T], ln_ext[l,t] <=  ln_s_nom[l]
        [l=1:LN_ext,t=1:T], ln_ext[l,t] >= -ln_s_nom[l]

        [l=1:LK_ext,t=1:T], lk_ext[l,t] >= lk_p_nom[l].*links_ext[l, :p_min_pu]
        [l=1:LK_ext,t=1:T], lk_ext[l,t] <= lk_p_nom[l].*links_ext[l, :p_max_pu]
end)


to_symbol(str) = Symbol(replace(str, " ", "_"))
#   nodal balance
@constraint(m, balance[n=1:N, t=1:T], (
      sum(g_fix[idx_by(generators_fix, :bus, [reverse_busidx[n]]), t])
    + sum(g_ext[idx_by(generators_ext, :bus, [reverse_busidx[n]]), t])
    # + sum(gcom[idx_by(generators_com, :bus, [reverse_busidx[1]]), t])

    - network.loads_t["p"][t,to_symbol(reverse_busidx[n])]

    + sum(ln_fix[ idx_by(lines_fix, :bus1, [reverse_busidx[n]]) ,t])
    - sum(ln_fix[ idx_by(lines_fix, :bus0, [reverse_busidx[n]]) ,t])
    + sum(ln_ext[ idx_by(lines_ext, :bus1, [reverse_busidx[n]]) ,t])
    - sum(ln_ext[ idx_by(lines_ext, :bus0, [reverse_busidx[n]]) ,t])

    + sum(lk_fix[ idx_by(links_fix, :bus1, [reverse_busidx[n]]) ,t])
    - sum(lk_fix[ idx_by(links_fix, :bus0, [reverse_busidx[n]]) ,t])
    + sum(lk_ext[ idx_by(links_ext, :bus1, [reverse_busidx[n]]) ,t])
    - sum(lk_ext[ idx_by(links_ext, :bus0, [reverse_busidx[n]]) ,t])
      == 0 ))



# Kirchhoff Voltage Law
for (branch, var, attribute) in [(lines, ln, :x), (links, lk, :r)]
    g = DiGraph(length(busidx))
    for l = eachrow(branch)
      add_edge!(g, busidx[l[:bus0]], busidx[l[:bus1]])
    end
    cycles = [i for i in simplecycles(g) if length(i)>2]
    if attribute==:x
        @constraint(m, line_cycle_constraint[c=1:length(cycles), t=1:T] ,
                sum(dot(lines[cycles[c], :x]/380. , ln[cycles[c],t])) == 0)
    end
    if attribute==:r
        @constraint(m, link_cycle_constraint[c=1:length(cycles), t=1:T] ,
                sum(dot(links[cycles[c], :r]/380. , lk[cycles[c],t])) == 0)
    end
end
# Might be nessecary to loop over all subgraphs as
# for (sn, sub) in enumerate(weakly_connected_components(g))
#     # g_sub = induced_subgraph(g, sub)[1]


@objective(m, Min, sum(dot(generators[:marginal_cost], gn[:,t]) for t=1:T))


solve(m)
@show getvalue(ln)
