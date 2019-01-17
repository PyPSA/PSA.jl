
using PSA, JuMP, Gurobi, AxisArrays, Revise
# include("/home/vres/data/lisa/jl/PSA.jl/src/auxilliaries.jl")
# n = PSA.import_nc("/home/lisa/Documents/Master-Arbeit/Code/test/three_nodes.nc")
#include("/home/lisa/vres/jl/PSA.jl/src/auxilliaries.jl")

# nname = "elmod_8760h"
# f_o = "1e3"
# n = PSA.import_nc("/home/vres/data/pypsa_models/elmod/$(nname).nc")
# PSA.calculate_dependent_values!(n, f_o)
# PSA.set_snapshots!(n, n.snapshots[1:10])

n = PSA.import_nc("/home/fabian/playground/pathway/elec_s_37_lv2.0_2020-2030_8760H.nc")
n.generators = PSA.assign(n.generators, 1e6, "maintenance_cost")
n.generators[:, "p_nom"] .= 0
n.generators[:, "p_nom_opt"] .= 0



# solver = GurobiSolver(Crossover=0)
# m = PSA.lopf_pathway(n, solver, investment_period=n.snapshots, 
#                      invest_at_first_sn=true)