using DataFrames, JuMP

import JuPSA
reload("JuPSA")
#
#
network = JuPSA.import_network("mininetwork")
network.loads_t["p"] = DataFrame(My_bus_0 = 0, My_bus_1 = 0, My_bus_2 = 100)
network.lines[3, :bus0], network.lines[3, :bus1] = network.lines[3, :bus1], network.lines[3, :bus0]
#
# network = JuPSA.import_network("/home/fabian/vres/py/pypsa/examples/ac-dc-meshed/ac-dc-data/")
# for name = [JuPSA.to_symbol(x) for x in network.buses[:name]]
#  if !in(name, names(network.loads_t["p"]))
#   network.loads_t["p"][name] = 0
#  end
# end


network = JuPSA.import_network("pre8-37")
network.loads_t["p"][:LT0_0] =0
# network.generators = network.generators[network.generators[:carrier].=="OCGT",:]
network.generators[:p_nom] = network.generators[:p_nom] *100
JuPSA.set_snapshots(network, 1:10)


JuPSA.calculate_dependent_values(network)
JuPSA.lopf(network)
println(network.lines_t["p0"])
