
using PSA, JuMP, Gurobi, AxisArrays, Revise
# include("/home/vres/data/lisa/jl/PSA.jl/src/auxilliaries.jl")
# n = PSA.import_nc("/home/lisa/Documents/Master-Arbeit/Code/test/three_nodes.nc")
#include("/home/lisa/vres/jl/PSA.jl/src/auxilliaries.jl")

nname = "elmod_8760h"
f_o = "1e3"

n = PSA.import_nc("/home/vres/data/pypsa_models/elmod/$(nname).nc")
PSA.calculate_dependent_values!(n, f_o)

# limit snapshots for testing
PSA.set_snapshots!(n, n.snapshots[1:2])

# lines which needs to be extended
extend_line = [676, 619, 264]
for line in extend_line
    n.lines[line,"capital_cost"] = 1e5/float(f_o)
    n.lines[line, "s_nom_extendable"] = 1
    n.lines[line,"s_nom_min"] = n.lines[line,"s_nom"]
end

# make lines extendable for testing
#= n.lines[:,"capital_cost"] = 1e10
n.lines[:, "s_nom_extendable"] = true
n.lines[:,"s_nom_min"] = n.lines[:,"s_nom"] =#

# set fixed p_nom
PSA.set_nuclear_p_nom!(n)
logfile = "/home/lisa/gurobi_log/$(nname)_$(f_o).txt"
solver = GurobiSolver(Crossover=0, LogFile=logfile)

# set CO2 limit
#n.global_constraints = AxisArray(fill(NaN, 6)', ["CO2Limit"], n.global_constraints.axes[2])
#n.global_constraints["CO2Limit", "constant"] = 1e5

PSA.lopf_pathway(n, solver, investment_period="y")