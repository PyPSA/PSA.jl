using PSA
using Test

const testdir = dirname("test/")

tests = [
    "PSA",
    "utils",
    "build_lopf",
    "build_block_lopf",
    "run_lopf",
    "run_iterative_lopf",
    "run_benders_lopf",
    "run_lazybenders_lopf"
]

@testset "PSA" begin
    for t in tests
        tp = joinpath(testdir, "$(t).jl")
        include(tp)
    end
end