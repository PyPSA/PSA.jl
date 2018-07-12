using PSA
using Test

const testdir = dirname("test/")

tests = [
    "auxiliaries",
    "lopf",
    "PSA",
    "netcdf_io"
]

@testset "PSA" begin
    for t in tests
        tp = joinpath(testdir, "$(t).jl")
        include(tp)
    end
end