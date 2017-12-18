import PSA
# reload("PSA")

network = PSA.import_network("/home/fabian/Desktop/PCA/intermediates/03-clusters/pre6-37/")

# network.loads_t["p"][:LT0_0] = 0
PSA.set_snapshots!(network, 1:100)
m = PSA.lopf(network; Crossover=1)
# println(network.storage_units_t["p"])
