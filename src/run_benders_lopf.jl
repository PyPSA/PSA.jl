# TODO: develop in IPYNB

using JuMP
#using CPLEX
using MathProgBase

include("utils.jl")

# set parameters
iteration = 1
tolerance = 1e-4

# set start solution

# run benders decomposition
function run_benders_lopf(network, solver)
    model_master = build_lopf(network, solver, benders="master")
    model_slave = build_lopf(network, solver, benders="slave")

    while(true) # set termination condition later

        status_master = solve(model_master)

        # cases of master problem
        if status_master == :Infeasible
            println("The problem is infeasible.")
            break
        elseif status_master == :Unbounded
            # TODO:
            fmCurrent = nothing
            xCurrent = nothing
        elseif status_master == :Optimal
            fmCurrent = getobjectivevalue(model_master)
            xCurrent = getvalue(investment_variables)
        else
            error("Weird status of master problem!")
        end

        # if iteration == 1 really needed?, maybe need a good starting solution..

        



    end
end

# write results