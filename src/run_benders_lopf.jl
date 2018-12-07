using JuMP
using MathProgBase

# TODO: neglects storage units, storages, and links (!) at the moment!
function run_benders_lopf(network, solver; 
                            formulation::String = "angles_linear",
                            investment_type::String = "continuous",
                            update_x::Bool = false, 
                            split::Bool = false,
                            individualcuts::Bool = false,
                            tolerance::Float64 = 1e-1,
                            mip_gap::Float64 = 1e-11,
                            bigM::Float64 = 1e12,
                            max_iterations::Int64 =  1000,
                            relax::Bool = false,
                            relax_threshold::Float64 = 1e5,
                            round::Bool = false,
                            round_threshold::Float64 = 0.5,
                            inexact_iterations::Int64 = 0,
                            inexact_mip_gap::Float64 = 1e-5,
                            inexact_threshold::Float64 = 1e5,
                            remove_inactive::Bool = false,
                            remove_inactive_threshold::Float64 = 1e2)

    calculate_dependent_values!(network)
    T = nrow(network.snapshots)
    ext_gens_b = ext_gens_b = convert(BitArray, network.generators[:p_nom_extendable])
    ext_lines_b = convert(BitArray, network.lines[:s_nom_extendable])
    N_ext_G = sum(ext_gens_b)
    N_ext_LN = sum(ext_lines_b)
    individualcuts ? N_groups = T : N_groups = 1

    coupled_slave_constrs = [
        :lower_gen_limit,
        :upper_gen_limit,
        :lower_line_limit,
        :upper_line_limit
        # TODO: add links
        ]
    
    model_master = build_lopf(network, solver,
        benders="master",
        investment_type=investment_type,
        N_groups=N_groups
    )

    models_slave = JuMP.Model[]
    if !split
        push!(models_slave, 
            build_lopf(network, solver, 
                benders="slave",
                formulation=formulation,
            )
        )
    else
        for t=1:T
            push!(models_slave, 
                build_lopf(network, solver, 
                    benders="slave",
                    formulation=formulation,
                    snapshot_number=t
                )
            )
        end
    end
    N_sub = length(models_slave)
        
    uncoupled_slave_constrs = setdiff(getconstraints(first(models_slave)), coupled_slave_constrs)
    vars_master = getvariables(model_master)
    vars_slave = getvariables(first(models_slave))
    
    p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=ext_gens_b)
    p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=ext_gens_b)
            
    iteration = 1
    inexact_iterations > 0 ? initial_mip_gap=inexact_mip_gap : initial_mip_gap = mip_gap
    push!(solver.options, (:MIPGap, initial_mip_gap))

    ubs = Float64[]
    lbs = Float64[]
    mem_master = Float64[] # MBytes

    println("\nITER\tGAP")

    while(iteration <= max_iterations)

        push!(mem_master,Base.summarysize(model_master)/1e6)
        
        status_master = solve(model_master, relaxation=relax);

        # cases of master problem
        if status_master == :Infeasible

            # println("The problem is infeasible.")
            break

        elseif status_master == :Unbounded

            objective_master_current = -bigM
            G_p_nom_current = network.generators[ext_gens_b,:][:p_nom]
            LN_s_nom_current = network.lines[ext_lines_b,:][:s_nom]
            
            # objective_master_current = bigM
            # G_p_nom_current = network.generators[ext_gens_b,:][:p_nom_max]
            # LN_s_nom_current = network.lines[ext_lines_b,:][:s_nom_max]
            
            LN_inv_current = zeros(nrow(network.lines[ext_lines_b,:]))
            setvalue(model_master[:ALPHA], -bigM)

        elseif status_master == :Optimal

            objective_master_current = getobjectivevalue(model_master)
            G_p_nom_current = getvalue(model_master[:G_p_nom])
            LN_s_nom_current = getvalue(model_master[:LN_s_nom])
            LN_inv_current = getvalue(model_master[:LN_inv])

            # inactive cut removal
            # TODO: can be supported with MathOptInterface

        else
            error("Odd status of master problem: $status_master")
        end

        if round && relax
            for l=1:N_ext_LN
                LN_s_nom_current[l] = (1.0+ceil(LN_inv_current[l])/network.lines[ext_lines_b,:num_parallel][l]) * network.lines[ext_lines_b,:s_nom][l]
            end
        end

        # # for debugging
        # println("Status of the master problem is ", status_master, 
        # "\nwith objective_master_current = ", objective_master_current, 
        # "\nG_p_nom_current = ", G_p_nom_current,
        # "\nLN_s_nom_current = ", LN_s_nom_current,
        # "\nalpha = ", getvalue(model_master[:ALPHA]))
        # investment_type=="integer" ? @show(LN_inv_current) : nothing
     
        # adapt RHS of slave model with solution from master problem
        model_slave = models_slave[1]
        for t=1:T

            split ? model_slave = models_slave[t] : nothing

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

        # if update_x
        #     network.lines[ext_lines_b,:][:x_pu] .= network.lines[ext_lines_b,:][:x_pu] .* network.lines[ext_lines_b,:][:s_nom] ./ LN_s_nom_current
        #     model_slave = build_lopf(network, solver, 
        #         benders="slave",
        #         formulation=formulation
        #     )
        # end

        statuses_slave = Symbol[]
        for i=1:N_sub
            push!(statuses_slave, solve(models_slave[i])) 
        end
        status_slave = minimum(statuses_slave)
        
        # # for debugging
        # if iteration == 1
        #     for i=1:length(models_slave)
        #         println(models_slave[i])
        #         @show(getvalue(models_slave[i][:G_ext]))
        #         @show(getvalue(models_slave[i][:G_fix]))
        #         @show(getobjectivevalue(models_slave[i]))
        #     end
        # end

        objective_slave_current = sum(getobjectivevalue(models_slave[i]) for i=1:N_sub)
        duals_lower_gen_limit = getduals(models_slave, :lower_gen_limit)
        duals_upper_gen_limit = getduals(models_slave, :upper_gen_limit)
        duals_lower_line_limit = getduals(models_slave, :lower_line_limit)
        duals_upper_line_limit = getduals(models_slave, :upper_line_limit)

        # # for debugging
        # println("Status of the slaveproblem is ", status_slave, 
        # "\nwith GAP ", objective_slave_current - getvalue(model_master[:ALPHA]),
        # "\nobjective_slave_current = ", objective_slave_current, 
        # "\nALPHA = ", getvalue(model_master[:ALPHA]))

        ALPHA_current = sum(getvalue(model_master[:ALPHA][g]) for g=1:N_groups)

        if objective_slave_current - ALPHA_current <= relax_threshold && relax == true
            println("Drop relaxation!")
            relax = false
        end

        if objective_slave_current - ALPHA_current <= inexact_threshold || inexact_iterations == iteration
            println("Drop inexact solutions!")
            println(mip_gap)
            push!(solver.options, (:MIPGap, mip_gap))
        end

        # cases of slave problem
        if (status_slave == :Optimal && 
            objective_slave_current - ALPHA_current <= tolerance)

            println("Optimal solution of the original problem found")
            println("The optimal objective value is ", objective_master_current)
            break

        end

        if !individualcuts
            Trange = [1:T] 
        else 
            Trange = 1:T
        end

        for Tr = Trange

            N_sub == 1 ? cnt = 1 : cnt = Tr
            individualcuts ? cntr = Tr : cntr = 1

            if (status_slave == :Optimal &&
                objective_slave_current - ALPHA_current > tolerance)

                # # for debugging
                # @show(get_benderscut_constant(models_slave[cnt],uncoupled_slave_constrs))
    
                # println("\nThere is a suboptimal vertex, add the corresponding constraint")
                @constraint(model_master, model_master[:ALPHA][cntr] >=
                
                    sum(  duals_lower_gen_limit[t,gr] * p_min_pu(t,gr) * model_master[:G_p_nom][gr]
                        + duals_upper_gen_limit[t,gr] * p_max_pu(t,gr) * model_master[:G_p_nom][gr]
                    for t=Tr for gr=1:N_ext_G)
                    + 
                    sum(  duals_lower_line_limit[t,l] * (-1) * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                        + duals_upper_line_limit[t,l] * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                    for t=Tr for l=1:N_ext_LN)
    
                    + get_benderscut_constant(models_slave[cnt],uncoupled_slave_constrs)
                
                )
            end
            
            if status_slave == :Infeasible
    
                # println("\nThere is an  extreme ray, adding the corresponding constraint")
                
                @constraint(model_master, 0 >=
    
                    sum(  duals_lower_gen_limit[t,gr] * p_min_pu(t,gr) * model_master[:G_p_nom][gr]
                        + duals_upper_gen_limit[t,gr] * p_max_pu(t,gr) * model_master[:G_p_nom][gr]
                    for t=Tr for gr=1:N_ext_G)
                    + 
                    sum(  duals_lower_line_limit[t,l] * (-1) * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                        + duals_upper_line_limit[t,l] * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                    for t=Tr for l=1:N_ext_LN)
    
                    + get_benderscut_constant(models_slave[cnt],uncoupled_slave_constrs)
                )
            end
            
        end

        push!(lbs, ALPHA_current)
        push!(ubs, objective_slave_current)
        println("$iteration\t$(objective_slave_current - ALPHA_current)")

        iteration += 1

    end

    if iteration <= max_iterations
        # converged
        write_optimalsolution(network, model_master; sm=models_slave, joint=false)
    else
        println("Hit the maximum number of iterations. No solution provided.\n
        Try choosing max_iterations>>$(max_iterations)!")
    end

    return (model_master, models_slave, lbs, ubs, mem_master);

end