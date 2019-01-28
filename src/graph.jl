using LightGraphs
using PyCall
networkx = pyimport("networkx" )

function to_graph(n)
    busidx = idx(n.buses)
    g = DiGraph(size(n.buses)[1])
    for l=1:size(n.lines)[1]
        add_edge!(g, busidx[n.lines[l,"bus0"]], busidx[n.lines[l,"bus1"]] )
    end
    for l=1:size(n.links)[1]
        add_edge!(g, busidx[n.links[l,"bus0"]], busidx[n.links[l,"bus1"]] )
    end
    return g
end

function to_simple_graph(n, branches=nothing)
    if branches == nothing
        @show(branches)
        branches = cat(n.lines[:, ["bus0", "bus1"]],
                        n.links[:, ["bus0", "bus1"]], dims=1)
    end
    busidx = idx(n.buses)
    g = SimpleGraph(size(n.buses)[1])
    for b=1:size(branches)[1]
        add_edge!(g, busidx[branches[b,"bus0"]], busidx[branches[b,"bus1"]] )
    end
    return g
end

function incidence_matrix(n)
    busidx = idx(n.buses)
    lines = n.lines
    K = zeros(size(n.buses)[1],size(lines)[1])
    for l in 1:size(K)[2]
        K[busidx[lines[l,"bus0"]],l] = 1
        K[busidx[lines[l,"bus1"]],l] = -1
    end
    return K
end

function laplace_matrix(n)
    K = incidence_matrix(n)
    return K*K'
end

function ptdf_matrix(n)
    K = incidence_matrix(n)
    H = K' * pinv(K*K')
    return H .- H[:,1]
end


function get_cycles(n, branches = nothing)
    branches == nothing ?  branches = n.lines : nothing
    LightGraphs.cycle_basis(to_simple_graph(n, branches))
end


# #returns connected buses in for each cycle in network
# function get_cycles(n, branches=nothing)
#     branches == nothing ?  branches = n.lines : nothing

#     busidx = idx(n.buses)
#     g = networkx[:Graph]()
#     g[:add_nodes_from](1:size(n.buses)[1])
#     g[:add_edges_from]([(busidx[branches[l,"bus0"]], busidx[branches[l,"bus1"]]) 
#                         for l=1:size(branches)[1]])
#     C = networkx[:cycle_basis](g)

#     #deal with multigraph:
#     busidx = idx(n.buses)
#     line_buses = sort(string.(n.lines[:, ["bus0", "bus1"]]), dims=2)
#     sorted_pairs = zip(line_buses[:, 1], line_buses[:, 2]) |> collect
#     for p in unique(sorted_pairs)
#         lines_i = findall((in)([p]), sorted_pairs)
#         if length(lines_i)>1
#             push!(C, [busidx[p[1]], busidx[p[2]]])
#         end
#     end
#     C
# end




#returns connected lines and their direction for each cycle
function get_directed_cycles(n, branches=nothing)

    branches == nothing ?  branches = n.lines : nothing

    L, = size(branches)
    rev_busidx = rev_idx(n.buses)
    cyc_bus = get_cycles(n, branches)

    cyc_bra = Array{Int64,1}[]
    dirs = Array{Float64,1}[]

    nor_order =  zip(string.(branches[:, "bus0"]), 
                     string.(branches[:, "bus1"])) |> collect
    rev_order =  zip(string.(branches[:, "bus1"]), 
                     string.(branches[:, "bus0"])) |> collect
    all_orders = [nor_order; rev_order]
    all_dirs = [fill(1., L); fill(-1., L)]


    function try_both_dirs(bus_tuple)
        if bus_tuple ∈ nor_order
            findall(in([bus_tuple]), nor_order)[1], 1. 
        else 
            findall(in([bus_tuple]), rev_order)[1], -1. 
        end
    end

    for (i_c, c) = enumerate(cyc_bus)
        if length(c) == 2
            bus0, bus1 = c[1], c[2]            
            bus0, bus1 = rev_busidx[bus0], rev_busidx[bus1]
            
            ix = findin(all_orders, [(bus0, bus1)])
            bs, ds = ix .% L, all_dirs[ix]
            push!(cyc_bra, bs)
            push!(dirs, ds)
                        
        else
            push!(cyc_bra, Int64[])
            push!(dirs, Float64[])
            for i_b = 1:length(c)

                bus0, bus1 = c[i_b], c[(i_b)%length(c) + 1]
                
                bus0, bus1 = rev_busidx[bus0], rev_busidx[bus1]
                @show(bus0, bus1)
                b, d = try_both_dirs((bus0, bus1))

                push!(cyc_bra[i_c], b)
                push!(dirs[i_c], d)
            end
        end
    end
    cyc_bra, dirs
end


# ------ fix_p_nom for nuclear plants -------------------
function set_nuclear_p_nom!(n)
    fix_p_nom = fill(NaN, (length(n.snapshots), length(n.generators[:, "p_nom_extendable"])))
    ext_gens_b = BitArray(n.generators[:, "p_nom_extendable"])

    nan_frame = (AxisArray(fix_p_nom, Axis{:time}(n.snapshots), 
                                   Axis{:col}(n.generators.axes[1][ext_gens_b])))
    replace_attribute!(n, :generators_t, "p_nom_opt", nan_frame)
    # add values 
    for i in 1:length(n.snapshots)
        if Dates.year(n.snapshots[i]) >= 2020
            n.generators_t.p_nom_opt[DateTime(n.snapshots[i]), "Nuclear 261"] = 0
        end
        if Dates.year(n.snapshots[i]) >= 2022
            n.generators_t.p_nom_opt[DateTime(n.snapshots[i]), "Nuclear 429"] = 0
            n.generators_t.p_nom_opt[DateTime(n.snapshots[i]), "Nuclear 217"] = 0
            n.generators_t.p_nom_opt[DateTime(n.snapshots[i]), "Nuclear 294"] = 0
        end
        if Dates.year(n.snapshots[i]) >= 2023
            n.generators_t.p_nom_opt[DateTime(n.snapshots[i]), "Nuclear 317"] = 0
            n.generators_t.p_nom_opt[DateTime(n.snapshots[i]), "Nuclear 392"] = 0
            n.generators_t.p_nom_opt[DateTime(n.snapshots[i]), "Nuclear 257"] = 0
        end
    end
    # return n
end
