using JuMP
using MathProgBase

"""Solves a linear optimal power flow"""
function run_lopf(network, solver; 
    rescaling::Bool=false,
    formulation::String="angles_linear",
    investment_type::String="continuous",
    blockmodel::Bool=false,
    decomposition::String="",
    n_sub::Int64=1)

    if blockmodel
        println("Build block JuMP model.")
        m = build_block_lopf(network, solver;
                rescaling=rescaling,
                formulation=formulation,
                investment_type=investment_type,
                decomposition=decomposition,
                n_sub=n_sub
            )
    else
        println("Build ordinary JuMP model.")
        m = build_lopf(network, solver; 
                rescaling=rescaling,
                formulation=formulation,
                investment_type=investment_type
            )
    end

    status = solve(m)

    if status==:Optimal
        write_optimalsolution(network, m)
    end

    return m

end