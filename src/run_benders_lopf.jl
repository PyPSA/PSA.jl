using JuMP
using MathOptInterface

# TODO: neglects storage units and  storages at the moment!

"""Solves a linear optimal power flow with Benders Decomposition."""
function run_benders_lopf(network, solver; 
    formulation::String = "angles_linear",
    investment_type::String = "continuous",
    rescaling::Bool = false,
    update_x::Bool = false, 
    split_subproblems::Bool = false, 
    individualcuts::Bool = false,
    tolerance::Float64 = 100.0,
    mip_gap::Float64 = 1e-8,
    bigM::Float64 = 1e12,
    max_iterations::Int64 =  1000,
    #relax::Bool = false, # not reimplemented in JuMP 0.19.0
    #relax_threshold::Float64 = 1e5,  # not reimplemented in JuMP 0.19.0
    round::Bool = false,
    round_threshold::Float64 = 0.5,
    inexact_iterations::Int64 = 0,
    inexact_mip_gap::Float64 = 1e-5,
    inexact_threshold::Float64 = 1e5,
    remove_inactive::Bool = false,
    remove_inactive_threshold::Float64 = 1e2,
    filter_duals::Bool = false
    )

    #################
    # sanity checks #
    #################

    # TODO: make a few checks on parameter selection!

    ###################
    # precalculations #
    ###################

    calculate_dependent_values!(network)
    T = nrow(network.snapshots)
    ext_gens_b = convert(BitArray, network.generators[:p_nom_extendable])
    ext_lines_b = convert(BitArray, network.lines[:s_nom_extendable])
    ext_links_b = convert(BitArray, network.links[:p_nom_extendable])
    N_ext_G = sum(ext_gens_b)
    N_ext_LN = sum(ext_lines_b)
    N_ext_LK = sum(ext_links_b)
    individualcuts ? N_cuts = T : N_cuts = 1
    
    p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=ext_gens_b)
    p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=ext_gens_b)
    filter_timedependent_extremes!(p_max_pu, 0.01)
    filter_timedependent_extremes!(p_min_pu, 0.01)
    candidates = line_extensions_candidates(network)

    bigm_upper = bigm(:flows_upper, network)
    bigm_lower = bigm(:flows_lower, network)

    rf = 1
    rf_dict = rescaling_factors(rescaling)

    coupled_slave_constrs = [
        :lower_bounds_G_ext,
        :upper_bounds_G_ext,
        :lower_bounds_LN_ext,
        :upper_bounds_LN_ext,
        ]

    if N_ext_LK > 0
        push!(coupled_slave_constrs, :lower_bounds_LK_ext)
        push!(coupled_slave_constrs, :upper_bounds_LK_ext)
    end

    if investment_type=="integer_bigm"
        push!(coupled_slave_constrs, :flows_upper)
        push!(coupled_slave_constrs, :flows_lower)
    end

    ###################################
    # build master and slave problems # 
    ###################################
    
    # build master model
    model_master = build_lopf(network, solver,
        benders="master",
        investment_type=investment_type,
        rescaling=rescaling,
        N_cuts=N_cuts
    )

    # build slave models
    models_slave = JuMP.Model[]
    if !split_subproblems
        push!(models_slave, 
            build_lopf(network, solver, 
                benders="slave",
                formulation=formulation,
                rescaling=rescaling,
            )
        )
    else
        for t=1:T
            push!(models_slave, 
                build_lopf(network, solver, 
                    benders="slave",
                    formulation=formulation,
                    rescaling=rescaling,
                    snapshot_number=t
                )
            )
        end
    end
    N_slaves = length(models_slave)
        
    uncoupled_slave_constrs = setdiff(getconstraints(first(models_slave)), coupled_slave_constrs)

    ############################
    # initialise log variables # 
    ############################

    iteration = 1
    inexact_iterations > 0 ? initial_mip_gap=inexact_mip_gap : initial_mip_gap = mip_gap
    inexact_iterations > 0 ? push!(solver.options, (:MIPGap, initial_mip_gap)) : nothing

    ubs = Float64[]
    lbs = Float64[]
    memory_master = Float64[] # MBytes

    println("\nITER\tSLAVE_OBJ_CURR\tALPHA_CURR\tGAP")

    ####################################
    # enter benders decomposition loop # 
    ####################################

    while(iteration <= max_iterations)

        push!(memory_master,Base.summarysize(model_master)/1e6)
       
        JuMP.optimize!(model_master)#, relaxation=relax);
        status_master = JuMP.termination_status(model_master)

        # cases of master problem
        if status_master == MOI.INFEASIBLE

            # print iis if infeasible
            if solver.constructor == Gurobi.Optimizer
                println("ERROR: Master problem is infeasible. The IIS is:")
                println(get_iis(model_master))
            end
            
            break

        elseif status_master == :Unbounded

            objective_master_current = -bigM
            G_p_nom_current = network.generators[ext_gens_b,:][:p_nom]
            LN_s_nom_current = network.lines[ext_lines_b,:][:s_nom]
            N_ext_LK>0 ? LK_p_nom_current = network.links[ext_links_b,:][:p_nom] : nothing
            investment_type=="integer_bigm" ? LN_opt_current = nothing : LN_inv_current = zeros(nrow(network.lines[ext_lines_b,:]))
            setvalue(model_master[:ALPHA], -bigM)

        elseif status_master == MOI.OPTIMAL

            objective_master_current = JuMP.objective_value(model_master)
            G_p_nom_current = JuMP.value.(model_master[:G_p_nom])
            LN_s_nom_current = JuMP.value.(model_master[:LN_s_nom])
            N_ext_LK>0 ? LK_p_nom_current = JuMP.value.(model_master[:LK_p_nom]) : nothing
            investment_type=="integer_bigm" ? LN_opt_current = JuMP.value.(model_master[:LN_opt]) : LN_inv_current = JuMP.value.(model_master[:LN_inv])

        else
            error("Odd status of master problem: $status_master")
        end

        if round && relax
            for l=1:N_ext_LN
                LN_s_nom_current[l] = (1.0+ceil(LN_inv_current[l])/network.lines[ext_lines_b,:num_parallel][l]) * network.lines[ext_lines_b,:s_nom][l]
            end
        end

        # if updatex create new subproblems
        if update_x
            network.lines[ext_lines_b,:][:x_pu] .= network.lines[ext_lines_b,:][:x_pu] .* network.lines[ext_lines_b,:][:s_nom] ./ LN_s_nom_current
            models_slave = JuMP.Model[]
            if !split_subproblems
                push!(models_slave, 
                    build_lopf(network, solver, 
                        benders="slave",
                        formulation=formulation,
                        rescaling=rescaling,
                    )
                )
            else
                for t=1:T
                    push!(models_slave, 
                        build_lopf(network, solver, 
                            benders="slave",
                            formulation=formulation,
                            rescaling=rescaling,
                            snapshot_number=t
                        )
                    )
                end
            end
        end
     
        # adapt RHS of slave model with solution from master problem
     
        model_slave = models_slave[1]
     
        for t=1:T

            split_subproblems ? model_slave = models_slave[t] : nothing

            # RHS generators

            rf = rf_dict[:bounds_G]

            for gr=1:N_ext_G

                JuMP.set_standard_form_rhs(
                    model_slave[:lower_bounds_G_ext][t,gr],
                    rf * G_p_nom_current[gr] * p_min_pu(t,gr)
                )

                JuMP.set_standard_form_rhs(
                    model_slave[:upper_bounds_G_ext][t,gr],
                    rf * G_p_nom_current[gr] * p_max_pu(t,gr)
                )

            end

            # RHS lines

            rf = rf_dict[:bounds_LN]

            for l=1:N_ext_LN

                JuMP.set_standard_form_rhs(
                    model_slave[:lower_bounds_LN_ext][t,l],
                    rf * (-1) * network.lines[ext_lines_b,:s_max_pu][l] * LN_s_nom_current[l]
                )
                
                JuMP.set_standard_form_rhs(
                    model_slave[:upper_bounds_LN_ext][t,l],
                    rf * network.lines[ext_lines_b,:s_max_pu][l] * LN_s_nom_current[l]
                )

            end

            # RHS links

            if N_ext_LK > 0

                rf = rf_dict[:bounds_LK]
    
                for l=1:N_ext_LK
    
                    JuMP.set_standard_form_rhs(
                        model_slave[:lower_bounds_LK_ext][t,l],
                        rf * LK_p_nom_current[l] * network.links[ext_links_b, :p_min_pu][l]
                    )
                    
                    JuMP.set_standard_form_rhs(
                        model_slave[:upper_bounds_LK_ext][t,l],
                        rf * LK_p_nom_current[l] * network.links[ext_links_b, :p_max_pu][l]
                    )
    
                end

            end

            # RHS flows
                    
            rf = rf_dict[:flows]

            for l=1:N_ext_LN
                
                if investment_type=="integer_bigm"
                   
                    for c in candidates[l]
                        
                        rhs = rf * ( LN_opt_current[l,c] - 1 ) * bigm_upper[l]
                        if (rhs < 1e-4 && rhs > -1e-4)
                            rhs = 0.0
                        end

                        JuMP.set_standard_form_rhs(model_slave[:flows_upper][l,c,t],rhs)
                        
                        rhs = rf * ( 1 - LN_opt_current[l,c] ) * bigm_lower[l]
                        if (rhs < 1e-4 && rhs > -1e-4)
                            rhs = 0.0
                        end

                        JuMP.set_standard_form_rhs(model_slave[:flows_lower][l,c,t], rhs)

                    end

                end

            end

        end

        # solve slave problems
        statuses_slave = MOI.TerminationStatusCode[]
        for i=1:N_slaves
            JuMP.optimize!(models_slave[i])
            push!(statuses_slave, JuMP.termination_status(models_slave[i])) 
        end
        status_slave = minimum(statuses_slave)

        # print iis if infeasible
        for i in length(statuses_slave)
            if statuses_slave[i] == MOI.INFEASIBLE && solver.constructor == Gurobi.Optimizer
                println("WARNING: Subproblem $i is infeasible. The IIS is:")
                println(get_iis(models_slave[i]))
            end
        end

        # get results of slave problems
        objective_slave_current = sum(JuMP.objective_value(models_slave[i]) for i=1:N_slaves)
        duals_lower_bounds_G_ext = duals(models_slave, :lower_bounds_G_ext, filter_b=filter_duals)
        duals_upper_bounds_G_ext = duals(models_slave, :upper_bounds_G_ext, filter_b=filter_duals)
        duals_lower_bounds_LN_ext = duals(models_slave, :lower_bounds_LN_ext, filter_b=filter_duals)
        duals_upper_bounds_LN_ext = duals(models_slave, :upper_bounds_LN_ext, filter_b=filter_duals)

        if N_ext_LK > 0
            duals_lower_bounds_LK_ext = duals(models_slave, :lower_bounds_LK_ext, filter_b=filter_duals)
            duals_upper_bounds_LK_ext = duals(models_slave, :upper_bounds_LK_ext, filter_b=filter_duals)
        end
        
        if investment_type=="integer_bigm"
            duals_flows_upper = duals_flows(models_slave, :flows_upper, filter_b=filter_duals)
            duals_flows_lower = duals_flows(models_slave, :flows_lower, filter_b=filter_duals)
        end

        ALPHA_current = sum(JuMP.value(model_master[:ALPHA][g]) for g=1:N_cuts)

        # relaxation handling
        # if objective_slave_current - ALPHA_current <= relax_threshold && relax == true
        #     println("Drop relaxation!")
        #     relax = false
        # end

        # inexactness handling
        if inexact_iterations!=0 && (objective_slave_current - ALPHA_current <= inexact_threshold || inexact_iterations == iteration)
            println("Drop inexact solutions!")
            println(mip_gap)
            push!(solver.options, (:MIPGap, mip_gap))
        end

        # go through cases of subproblems
        if (status_slave == MOI.OPTIMAL && 
            abs(objective_slave_current - ALPHA_current) <= tolerance)

            push!(lbs, ALPHA_current)
            push!(ubs, objective_slave_current)
            println("$iteration\t$objective_slave_current\t$ALPHA_current\t$(objective_slave_current - ALPHA_current)")

            println("Optimal solution of the original problem found")
            println("The optimal objective value is ", objective_master_current)
            break

        end

        if !individualcuts
            T_range = [1:T] 
        else 
            T_range = 1:T
        end

        for T_range_curr = T_range

            N_slaves == 1 ? slave_id = 1 : slave_id = T_range_curr
            individualcuts ? cut_id = T_range_curr : cut_id = 1

            # calculate coupled cut components

            cut_G_lower = sum(  
                    duals_lower_bounds_G_ext[t,gr] * rf_dict[:bounds_G] * p_min_pu(t,gr) * model_master[:G_p_nom][gr]
                for t=T_range_curr for gr=1:N_ext_G)
            
            cut_G_upper = sum(  
                    duals_upper_bounds_G_ext[t,gr] * rf_dict[:bounds_G] * p_max_pu(t,gr) * model_master[:G_p_nom][gr]
                for t=T_range_curr for gr=1:N_ext_G)

            cut_LN_lower = sum(  
                    duals_lower_bounds_LN_ext[t,l] * rf_dict[:bounds_LN] * (-1) * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                for t=T_range_curr for l=1:N_ext_LN)

            cut_LN_upper = sum(  
                    duals_upper_bounds_LN_ext[t,l] * rf_dict[:bounds_LN] * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                for t=T_range_curr for l=1:N_ext_LN)

            cut_LK_lower = 0
            cut_LK_upper = 0

            if N_ext_LK > 0
                cut_LK_lower = sum( 
                        duals_lower_bounds_LK_ext[t,l] * rf_dict[:bounds_LK] * network.links[ext_links_b,:p_min_pu][l] * model_master[:LK_p_nom][l]
                    for t=T_range_curr for l=1:N_ext_LK)
    
                cut_LK_upper = sum( 
                        duals_upper_bounds_LK_ext[t,l] * rf_dict[:bounds_LK] * network.links[ext_links_b,:p_max_pu][l] * model_master[:LK_p_nom][l]
                    for t=T_range_curr for l=1:N_ext_LK)
            end

            cut_G = cut_G_lower + cut_G_upper
            cut_LN = cut_LN_lower + cut_LN_upper
            cut_LK = cut_LK_lower + cut_LK_upper

            if investment_type=="integer_bigm"

                cut_flows_lower = sum(  
                        duals_flows_lower[slave_id][l,c+1,t] * rf_dict[:flows] *( 1 - model_master[:LN_opt][l,c] ) * bigm_lower[l] 
                    for t=T_range_curr for l=1:N_ext_LN for c in candidates[l])

                cut_flows_upper = sum(  
                        duals_flows_upper[slave_id][l,c+1,t] * rf_dict[:flows] * ( model_master[:LN_opt][l,c] - 1 ) * bigm_upper[l] 
                    for t=T_range_curr for l=1:N_ext_LN for c in candidates[l])

                cut_flows = cut_flows_lower + cut_flows_upper

            end

            # calculate uncoupled cut components
            cut_const = get_benderscut_constant(models_slave[slave_id],uncoupled_slave_constrs)

            if (status_slave == MOI.OPTIMAL &&
                abs(objective_slave_current - ALPHA_current) > tolerance)

                rf = rf_dict[:benderscut]

                if investment_type!="integer_bigm"

                    @constraint(model_master, rf * model_master[:ALPHA][cut_id] >= 
                        rf * ( cut_G + cut_LN + cut_LK + cut_const ) 
                    )

                else

                    @constraint(model_master, rf * model_master[:ALPHA][cut_id] >= 
                        rf * ( cut_G + cut_LN + cut_LK + cut_flows + cut_const ) 
                    )

                end
            end
            
            if status_slave == MOI.INFEASIBLE

                if investment_type!="integer_bigm"

                    @constraint(model_master, 0 >= 
                        rf * ( cut_G + cut_LN + cut_LK + cut_const ) 
                    )
                
                else
                    
                    @constraint(model_master, 0 >=
                        rf * ( cut_G + cut_LN + cut_LK + cut_flows + cut_const ) 
                    )

                end
                
            end

        end

        push!(lbs, ALPHA_current)
        push!(ubs, objective_slave_current)
        println("$iteration\t$objective_slave_current\t$ALPHA_current\t$(objective_slave_current - ALPHA_current)")

        iteration += 1

    end # benders decomposition loop

    ##################
    # output results #
    ##################

    if iteration <= max_iterations
        write_optimalsolution(network, model_master; sm=models_slave, joint=false)
    else
        println("WARNING: Hit the maximum number of iterations. No solution provided.\n
        Try setting max_iterations > $(max_iterations)!")
    end

    return (model_master, models_slave, lbs, ubs, memory_master);

end