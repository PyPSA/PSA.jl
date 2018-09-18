using JuMP
using CPLEX
using MathProgBase

include("utils.jl")

# set parameters
iteration = 1
tolerance = 1e-4

# set start solution

# set master problem formulation
function build_master_lopf(network, solver)

end

# set slave problem formulation
function build_slave_lopf(network, solver)

end

# run benders decomposition
function run_lopf(network, solver)
    model_master = build_master_lopf(network, solver)
    model_slave = build_slave_lopf(network, solver)

    while(true) # set termination condition later

        status_master = solve(model_master)

        # cases of master problem
        if status_master == :Infeasible

        elseif status_master == :Unbounded

        elseif status_master == :Optimal

        else
            error("Weird status of master problem!")
        end

        # if iteration == 1 really needed?, maybe need a good starting solution..

        



    end
end

# write results