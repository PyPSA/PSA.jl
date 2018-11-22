using JuMP
using MathProgBase

function run_lopf(network, solver; rescaling::Bool=false,formulation::String="angles_linear", objective::String="total", investment_type::String="continuous", blockmodel::Bool=false, decomposition::String="")

    if blockmodel
        println("Build block JuMP model.")
        m = build_block_lopf(network, solver; formulation=formulation, investment_type=investment_type,decomposition=decomposition)
    else
        println("Build ordinary JuMP model.")
        m = build_lopf(network, solver; rescaling=rescaling,formulation=formulation, objective=objective, investment_type=investment_type)
    end

    status = solve(m)

    if status==:Optimal
        println(typeof(network))
        write_optimalsolution(network, m)
    end

    return m

end