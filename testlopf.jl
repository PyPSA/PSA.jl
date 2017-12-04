import JuPSA
reload("JuPSA")

network = JuPSA.import_network("/home/fabian/vres/py/pypsa/examples/ac-dc-meshed/ac-dc-data/")
for name = [JuPSA.to_symbol(x) for x in network.buses[:name]]
    if !in(name, names(network.loads_t["p"]))
        network.loads_t["p"][name] = 0
    end
end

# network = JuPSA.import_network("/home/fabian/Desktop/PCA/intermediates/03-clusters/pre6-37")
# JuPSA.set_snapshots!(network, 1:10)

# 

# network.loads_t["p"][:LT0_0] = 0
m = JuPSA.lopf(network; Crossover=0)
println(network.storage_units_t["p"])
