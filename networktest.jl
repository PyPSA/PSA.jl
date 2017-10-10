using DataFrames, JuMP, LightGraphs, Gurobi, PyCall
using JuPSA

@pyimport networkx



# network = import_network("mininetwork")
# network.loads_t["p"] = DataFrame(My_bus_0 = 0, My_bus_1 = 0, My_bus_2 = 100)
# network.lines[3, :bus0], network.lines[3, :bus1] = network.lines[3, :bus1], network.lines[3, :bus0]

network = JuPSA.import_network("pre8-37")
network.loads_t["p"][:LT0_0] =0
network.generators = network.generators[network.generators[:carrier].=="OCGT",:]
network.generators[:p_nom] = network.generators[:p_nom] *100
network.loads_t["p"] = network.loads_t["p"][1:5,:]

JuPSA.calculate_dependent_values(network)
# network.lines[1, :extendable] = true
# network.lines[1, :s_nom_min] = 0
# network.lines[1, :s_nom_max] = 1000
# network.generators[1,:commitable] = true
network.generators[:min_up_time] = 0
# network.generators[1,:min_up_time] = 5
network.generators[:min_down_time] = 0
# network.generators[1,:min_down_time] = 5
network.generators[:initial_status] = true


m = Model(solver=GurobiSolver())

buses = network.buses
reverse_busidx = rev_idx(buses)
busidx = idx(buses)
N = nrow(network.buses)
T = nrow(network.loads_t["p"]) #normally snapshots


# 1. add all generators to the model
# 1.1 set different generator types
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
append_idx_col([generators_fix, generators_ext, generators_com, generators])

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
if length(network.lines)>0
    lines_fix = network.lines[!network.lines[:extendable],:]
    lines_ext = network.lines[network.lines[:extendable],:]
    lines = [lines_fix; lines_ext]
else
    lines_fix = network.lines; lines_ext = network.lines; lines = network.lines
end
append_idx_col([lines_fix, lines_ext, lines])

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
if length(network.links)>0
    links_fix = network.links[!network.links[:extendable],:]
    links_ext = network.links[network.links[:extendable],:]
    links = [links_fix; links_ext]
else
    links_fix = network.links; links_ext = network.links; links = network.links
end
append_idx_col([links_fix, links_ext, links])

# 3.2 iterator bounds
LK_fix = nrow(links_fix)
LK_ext = nrow(links_ext)

#  3.3 set link variables
@variables m begin
    -links_fix[l, :p_nom].*links_fix[l, :p_min_pu]  <=  lk_fix[l=1:LK_fix,t=1:T]        <= links_fix[l, :p_nom].*links_fix[l, :p_max_pu]
    lk_ext[l=1:LK_ext,t=1:T]
    links_ext[l, :p_nom_min]                <=  lk_p_nom[l=1:LK_ext]     <= links_ext[l, :p_nom_max]
end
lk = [lk_fix; lk_ext]

# 3.4 set constraints for extendable links
@constraints(m, begin
        [l=1:LK_ext,t=1:T], lk_ext[l,t] >= lk_p_nom[l].*links_ext[l, :p_min_pu]
        [l=1:LK_ext,t=1:T], lk_ext[l,t] <= lk_p_nom[l].*links_ext[l, :p_max_pu]
end)


# 4. define nodal balance constraint
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


# 5. set Kirchhoff Voltage Law constraint
# since cyclebasis is not yet supported in LightGraphs, use the pyimported
# netwrokx in order to define all cycles. The cycle_basis returns a list of
# cycles, indicating the connected buses. For each cycle the connecting branches
# and their directions (signs) have to be determined.

# Might be nessecary to loop over all subgraphs as
# for (sn, sub) in enumerate(weakly_connected_components(g))
#     # g_sub = induced_subgraph(g, sub)[1]


for (branch, var, attribute) in [(lines, ln, :x)]#, (links, lk, :r)]
    g = networkx[:Graph]()
    g[:add_nodes_from](busidx)
    g[:add_edges_from]([(busidx[l[:bus0]], busidx[l[:bus1]]) for l in eachrow(lines)])
    cycles = [cycle for cycle in networkx[:cycle_basis](g) if length(cycle)>2]
    if length(cycles)>0
        cycles_branch = Array{Int64,1}[]
        directions = Array{Float64,1}[]
        for cyc=1:length(cycles)
            push!(cycles_branch,Int64[])
            push!(directions,Float64[])
            for bus=1:length(cycles[cyc])
                bus0 = cycles[cyc][bus]
                if bus == length(cycles[cyc])
                    bus1 = cycles[cyc][1]
                else
                    bus1 = cycles[cyc][bus+1]
                end
                try
                    push!(cycles_branch[cyc], branch[((branch[:bus0].==reverse_busidx[bus0])
                                .&(branch[:bus1].==reverse_busidx[bus1])),:idx][1] )
                    push!(directions[cyc], 1.)
                catch y
                    if isa(y, BoundsError)
                        push!(cycles_branch[cyc], branch[((branch[:bus0].==reverse_busidx[bus1])
                                        .&(branch[:bus1].==reverse_busidx[bus0])),:idx][1] )
                        push!(directions[cyc], -1.)
                    else
                        return y
                    end
                end
            end
        end
        if attribute==:x
            @constraint(m, line_cycle_constraint[c=1:length(cycles_branch), t=1:T] ,
                    dot(directions[c] .* lines[cycles_branch[c], :x]/380. , ln[cycles_branch[c],t]) == 0)
        end
        if attribute==:r
            @constraint(m, link_cycle_constraint[c=1:length(cycles_branch), t=1:T] ,
                    dot(directions[c] .* links[cycles_branch[c], :r]/380. , lk[cycles_branch[c],t]) == 0)
        end
    end
end


# 6. set objective function
@objective(m, Min, sum(dot(generators[:marginal_cost], gn[:,t]) for t=1:T))


solve(m)

# 7. extract optimisation results
network.generators_t["p"] = names!(DataFrame(transpose(getvalue(gn))),
                    [to_symbol(n) for n=generators[:name]])
network.lines_t["p0"] = names!(DataFrame(transpose(getvalue(ln))),
                    [to_symbol(n) for n=lines[:name]])

network.links_t["p0"] = names!(DataFrame(transpose(getvalue(lk))),
                    [to_symbol(n) for n=links[:name]])
