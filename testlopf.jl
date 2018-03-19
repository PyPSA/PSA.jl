import PSA
using DataFrames
# reload("PSA")

# network = PSA.import_network("/home/fabian/pca/intermediates/03-clusters/pre6-37")
network = PSA.import_network("/home/fabian/playground/37");

PSA.set_snapshots!(network, 1:10)
m = PSA.lopf(network; Crossover=1)
# println(network.storage_units_t["p"])
println(network.generators_t["p"])