
# Cf. PyPSA example https://github.com/PyPSA/PyPSA/tree/master/examples/ac-dc-meshed

push!(LOAD_PATH, "/home/tom/fias/lib")

import PSA

using Clp

network = PSA.import_network("/home/tom/fias/lib/pypsa/examples/ac-dc-meshed/ac-dc-data/")

solver = ClpSolver()

m = PSA.lopf(network, solver)

print(m.objVal)

print(network.generators_t["p"])
