using DataFrames, JuMP
using NetCDF

import JuPSA
reload("JuPSA")
#
#
# network = JuPSA.import_network("mininetwork")
# network.loads_t["p"] = DataFrame(My_bus_0 = 0, My_bus_1 = 0, My_bus_2 = 100)
# network.lines[3, :bus0], network.lines[3, :bus1] = network.lines[3, :bus1], network.lines[3, :bus0]
#
# network = JuPSA.import_network("/home/fabian/vres/py/pypsa/examples/ac-dc-meshed/ac-dc-data/")
# for name = [JuPSA.to_symbol(x) for x in network.buses[:name]]
#  if !in(name, names(network.loads_t["p"]))
#   network.loads_t["p"][name] = 0
#  end
# end

network = JuPSA.import_network("pre8-37")


# NetCDF.create("network.nc")


network.loads_t["p"][:LT0_0] = 0
# JuPSA.calculate_dependent_values(network)
# JuPSA.set_snapshots(network, 1:10)
# #
# m = JuPSA.lopf(network)
# println(network.storage_units_t["p"])
