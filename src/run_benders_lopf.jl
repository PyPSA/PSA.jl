using JuMP
#using CPLEX
using MathProgBase

# set parameters
tolerance = 1e-12 # TODO: adapt to needs
bigM = 1e12 # TODO: high enough?
max_iterations = 30

# run benders decomposition
# TODO: neglects storage units, storages, and links (!) at the moment!
function run_benders_lopf(network, solver)

    calculate_dependent_values!(network)
    T = nrow(network.snapshots)
    ext_gens_b = ext_gens_b = convert(BitArray, network.generators[:p_nom_extendable])
    ext_lines_b = convert(BitArray, network.lines[:s_nom_extendable])
    N_ext_G = sum(ext_gens_b)
    N_ext_LN = sum(ext_lines_b)

    coupled_slave_constrs = [
        :lower_gen_limit,
        :upper_gen_limit,
        :lower_line_limit,
        :upper_line_limit
        # TODO: add other components later
        ]

        
    model_master = build_lopf(network, solver, benders="master")
    model_slave = build_lopf(network, solver, benders="slave")
        
    uncoupled_slave_constrs = setdiff(getconstraints(model_slave), coupled_slave_constrs)
    
    vars_master = getvariables(model_master)
    vars_slave = getvariables(model_slave)
    
    p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=ext_gens_b)
    p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=ext_gens_b)
            
    iteration = 1

    while(iteration <= max_iterations)
        
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
            error("Odd status of master problem: $status_master")
        end

        println("Status of the master problem is ", status_master, 
        "\nwith objective_master_current = ", objective_master_current, 
        "\nG_p_nom_current = ", G_p_nom_current,
        "\nLN_s_nom_current = ", LN_s_nom_current,
        "\nalpha = ", getvalue(model_master[:ALPHA]))
     
        # adapt RHS of slave model with solution from master problem
        for t=1:T

            for gr=1:N_ext_G
                JuMP.setRHS(
                    model_slave[:lower_gen_limit][t,gr],
                    G_p_nom_current[gr]*p_min_pu(t,gr)
                )
                JuMP.setRHS(
                    model_slave[:upper_gen_limit][t,gr],
                    G_p_nom_current[gr]*p_max_pu(t,gr)
                )
            end

            for l=1:N_ext_LN
                JuMP.setRHS(
                    model_slave[:lower_line_limit][t,l],
                    -network.lines[ext_lines_b,:s_max_pu][l]*LN_s_nom_current[l]
                )
                JuMP.setRHS(
                    model_slave[:upper_line_limit][t,l],
                    network.lines[ext_lines_b,:s_max_pu][l]*LN_s_nom_current[l]
                )
            end

        end

        status_slave = solve(model_slave)

        investment_cost = objective_master_current - getvalue(model_master[:ALPHA])
        objective_slave_current = investment_cost + getobjectivevalue(model_slave)
        duals_lower_gen_limit = getdual(model_slave[:lower_gen_limit])
        duals_upper_gen_limit = getdual(model_slave[:upper_gen_limit])
        duals_lower_line_limit = getdual(model_slave[:lower_line_limit])
        duals_upper_line_limit = getdual(model_slave[:upper_line_limit])

        println("Status of the slaveproblem is ", status_slave, 
        "\nwith GAP ", objective_slave_current - objective_master_current,
        "\nobjective_slave_current = ", objective_slave_current, 
        "\nobjective_master_current = ", objective_master_current)

        # cases of slave problem
        if (status_slave == :Optimal && 
            #objective_slave_current - objective_master_current <= tolerance)
            objective_slave_current / objective_master_current <= 1 + tolerance)

            println("Optimal solution of the original problem found")
            println("The optimal objective value is ", objective_master_current)
            break

        end
        
        if (status_slave == :Optimal &&
            objective_slave_current > objective_master_current)

            println("\nThere is a suboptimal vertex, add the corresponding constraint")
            
            @constraint(model_master, model_master[:ALPHA] >=

                sum(  duals_lower_gen_limit[t,gr] * (-1) * p_min_pu(t,gr) * model_master[:G_p_nom][gr]
                    + duals_upper_gen_limit[t,gr] * p_max_pu(t,gr) * model_master[:G_p_nom][gr]
                for t=1:T for gr=1:N_ext_G)
                + 
                sum(  duals_lower_line_limit[t,l] * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                    + duals_upper_line_limit[t,l] * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                for t=1:T for l=1:N_ext_LN)

                + get_benderscut_constant(model_slave,uncoupled_slave_constrs)
            
            )

        end
        
        if status_slave == :Infeasible

            println("\nThere is an  extreme ray, adding the corresponding constraint")
            
            @constraint(model_master, 0 >=

                sum(  duals_lower_gen_limit[t,gr] * (-1) * p_min_pu(t,gr) * model_master[:G_p_nom][gr]
                    + duals_upper_gen_limit[t,gr] * p_max_pu(t,gr) * model_master[:G_p_nom][gr]
                for t=1:T for gr=1:N_ext_G)
                + 
                sum(  duals_lower_line_limit[t,l] * lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                    + duals_upper_line_limit[t,l] * lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                for l=1:N_ext_LN)

                + get_benderscut_constant(model_slave,uncoupled_slave_constrs)
            
            )
        end
       
        iteration += 1

    end
end

# write results