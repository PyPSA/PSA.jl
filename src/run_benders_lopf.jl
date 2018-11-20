using JuMP
#using CPLEX
using MathProgBase

include("utils.jl")

# set parameters
iteration = 1
tolerance = 1e-4
bigM = 1e12 # TODO: high enough?

# run benders decomposition
# TODO: neglects storage units, storages, and links at the moment!
function run_benders_lopf(network, solver)

    T = nrow(snapshots)
    ext_gens_b = ext_gens_b = convert(BitArray, generators[:p_nom_extendable])
    ext_lines_b = convert(BitArray, network.lines[:s_nom_extendable])
    N_ext_G = sum(ext_gens_b)
    N_ext_LN = sum(ext_lines_b)

    constraint_update_names = [
        :lower_gen_limit,
        :upper_gen_limit,
        :lower_line_limit,
        :upper_line_limit
        # TODO: add other components later
        ]
        
    model_master = build_lopf(network, solver, benders="master")
    model_slave = build_lopf(network, solver, benders="slave")

    vars_master = getvariables(model_master)
    vars_slave = getvariables(model_slave)
            
    while(true) # set termination condition later
        
        println("\n-----------------------")
        println("Iteration number = ", iteration)
        println("-----------------------\n")
        status_master = solve(model_master)

        # cases of master problem
        if status_master == :Infeasible

            println("The problem is infeasible.")
            break

        elseif status_master == :Unbounded

            objective_master_current = bigM
            G_p_nom_current = network.generators[ext_gens_b,:][:p_nom_max]
            LN_s_nom_current = network.lines[ext_lines_b,:][:s_nom_max]

        elseif status_master == :Optimal

            objective_master_current = getobjectivevalue(model_master)
            G_p_nom_current = getvalue(model_master[:G_p_nom])
            LN_s_nom_current = getvalue(model_master[:LN_s_nom])

        else
            error("Odd status of master problem!")
        end

        println("Status of the master problem is ", status_master, 
        "\nwith objective_master_current = ", objective_master_current, 
        "\nG_p_nom_current = ", G_p_nom_current,
        "\nLN_s_nom_current = ", LN_s_nom_current,
        "\nalpha = ", getvalue(model[:ALPHA]))
     
        # adapt RHS of slave model with solution from master problem

        for t=1:T

            for gr=1:N_ext_G
                JuMP.setRHS(
                    model[:lower_gen_limit][t,gr],
                    G_p_nom_current[gr]*resc_factor*p_min_pu(t,gr)
                )
                JuMP.setRHS(
                    model[:upper_gen_limit][t,gr],
                    G_p_nom_current[gr]*resc_factor*p_max_pu(t,gr)
                )
            end

            for l=1:N_ext_LN
                JuMP.setRHS(
                    model[:lower_line_limit][t,l],
                    -resc_factor*lines[:s_max_pu][l]*LN_s_nom_current[l]
                )
                JuMP.setRHS(
                    model[:upper_line_limit][t,l],
                    resc_factor*lines[:s_max_pu][l]*LN_s_nom_current[l]
                )
            end

        end

        print("\nThe current slave-problem model is \n", model_slave)

        status_slave = solve(model_slave)

        objective_slave_current = getobjectivevalue(model_slave)
        duals_lower_gen_limit = getdual(model[:lower_gen_limit])
        duals_upper_gen_limit = getdual(model[:upper_gen_limit])
        duals_lower_line_limit = getdual(model[:lower_line_limit])
        duals_upper_line_limit = getdual(model[:upper_line_limit])

        println("Status of the slaveproblem is ", status_slave, 
        "\nwith GAP ", objective_slave_current - objective_master_current,
        "\nobjective_slave_current = ", objective_slave_current, 
        "\nobjective_master_current = ", objective_master_current)

        # cases of slave problem
        if (status_slave == :Optimal && 
            objective_slave_current - objective_master_current <= tolerance)

            println("Optimal solution of the original problem found")
            println("The optimal objective value is ", objective_master_current)
            break

        end
        
        if (status_slave == :Optimal &&
            objective_slave_current > objective_master_current)

            println("\nThere is a suboptimal vertex, add the corresponding constraint")
            # TODO: add cuts
            #@constraint(model_master,alpha>=(muCurrent.*supply)'*x + sum(lambdaCurrent.*demand))
        
        end
        
        if status_slave == :Infeasible

            println("\nThere is an  extreme ray, adding the corresponding constraint")
            # TODO: add cuts
            #@constraint(model_master,0>=(muCurrent.*supply)'*x + sum(lambdaCurrent.*demand))
        
        end
        
        iteration += 1

    end
end

# write results