# Imports
import Base.Test: @inferred
using BenchmarkTools
using LightGraphs
using PyCall
const networkx = PyNULL()
copy!(networkx, pyimport("networkx"))

# Generate random graph with networkx
Gx = networkx[:fast_gnp_random_graph](1000, 0.5)
nodes = Gx[:nodes]
edges = Gx[:edges];

# Get this graph into LightGraphs
G = Graph()
add_vertices!(G, length(nodes))
for e in edges add_edge!(G, e) end

# Analyse
@btime cycle_basis($G, nothing)

@code_warntype cycle_basis(G, nothing)
