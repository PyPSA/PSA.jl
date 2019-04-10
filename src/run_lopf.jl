using JuMP
using MathOptInterface
using PowerModels
pm = PowerModels

"""Solves a linear optimal power flow"""
function run_lopf(network, solver; 
    rescaling::Bool=false,
    formulation="angles_linear",#::Union{String, Type{GenericPowerModel{F}}}="angles_linear",
    investment_type::String="continuous",
    blockmodel::Bool=false,
    decomposition::String="",
    n_sub::Int64=1) #where F <: pm.AbstractPowerFormulation

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

    JuMP.optimize!(m)
    status = JuMP.termination_status(m)

    if status==MOI.OPTIMAL
        write_optimalsolution(network, m)
    end

    # print iis if infeasible
    if status == MOI.INFEASIBLE && solver.constructor == Gurobi.Optimizer
        println("WARNING: Problem is infeasible. The IIS is:")
        println(get_iis(m))
    end

    return m

end