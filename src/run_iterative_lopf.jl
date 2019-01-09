"""Solves an iterative linear optimal power flow"""
function run_iterative_lopf(network, solver, iterations; 
    rescaling::Bool=false,
    formulation::String="angles_linear",
    investment_type::String="continuous",
    post_discretization::Bool=false,
    discretization_thresholds::Array{Float64,1}=[0.2,0.3],
    blockmodel=false,
    decomposition=""
    )

    # precalculations and initialisation
    x_0 = deepcopy(network.lines[:x])
    s_nom_0 = deepcopy(network.lines[:s_nom])
    objectives = Float64[]
    capacities = Array{Float64,1}[]
    reactances = Array{Float64,1}[]
    m = nothing

    # rounding function
    function round_line_extension!()
        for l=1:nrow(network.lines)

            if network.lines[:s_nom_extendable][l]

                num_parallel_extension = (s_nom_opt_continuous[l] / s_nom_0[l] - 1) * num_parallel_0[l]
                if mod(num_parallel_extension,1) >= threshold 
                    num_parallel_extension = ceil(num_parallel_extension)
                else 
                    num_parallel_extension = floor(num_parallel_extension)
                end

                extension_factor = (num_parallel_extension / num_parallel_0[l]+1)
                network.lines[:x][l] = x_0[l] / extension_factor
                network.lines[:s_nom_opt][l] = s_nom_0[l] * extension_factor
                network.lines[:s_nom][l] = network.lines[:s_nom_opt][l]
                network.lines[:s_nom_extendable][l] = false   
            
            end
        end
    end


    # solve
    for k=1:iterations

        println("Run iteration $k")

        if decomposition=="benders" && blockmodel==false
            m, sm = run_lazybenders_lopf(network, solver,
                rescaling=rescaling,
                formulation=formulation,
                investment_type=investment_type,
                split_subproblems=true,
                individualcuts=true)
        else
            m = run_lopf(network, solver; 
                rescaling=rescaling,
                formulation=formulation,
                investment_type=investment_type,
                blockmodel=blockmodel,
                decomposition=decomposition
                )
        end


        push!(objectives, m.objVal)
        push!(capacities, network.lines[:s_nom_opt])
        push!(reactances, network.lines[:x])

        for l=1:nrow(network.lines)
            if network.lines[:s_nom_extendable][l]
                if network.lines[:s_nom_opt][l] == 0.0
                    # reactance cannot take infinity values, instead choose prohibitively high value!
                    network.lines[:x][l] = 10e7 
                else
                    network.lines[:x][l] = (x_0[l] * s_nom_0[l]) / network.lines[:s_nom_opt][l]
                end
            end
        end

    end

    # perform post discretization if selected
    if post_discretization

        # store the optimal solution of continuous optimisation
        s_nom_opt_continuous = deepcopy(network.lines[:s_nom_opt])
        s_nom_extendable_0 = deepcopy(network.lines[:s_nom_extendable])
        num_parallel_0 = deepcopy(network.lines[:num_parallel])
        m_opt = nothing
        threshold = discretization_thresholds[1]
        threshold_opt = threshold

        # iterate through all possible discretization thresholds
        println("LENGTH OF DISC THRESH IS $(length(discretization_thresholds))")
        if (length(discretization_thresholds) > 1)

            for threshold in discretization_thresholds

                println("\nSTART EVALUATING THRESHOLD $threshold\n")

                # run lopf with rounded line capacities 
                round_line_extension!()
                m_threshold = run_lopf(network, solver; rescaling=rescaling)

                # compare to best solution in loop; better gets model
                if (threshold == discretization_thresholds[1]) || (m_threshold.objVal < m_opt.objVal)
                    println("\nCURRENT THRESHOLD BETTER OR FIRST ITERATION -- UPDATING threshold_opt to $threshold\n")
                    m_opt = m_threshold
                    threshold_opt = threshold
                end

                # reset network extendable circuits
                network.lines[:s_nom_extendable] = deepcopy(s_nom_extendable_0)

            end
        end

        # run with optimal threshold choice
        println("\nRUNNING AGAIN WITH OPTIMAL THRESHOLD CHOICE $threshold_opt\n")

        round_line_extension!()
        s_nom_opt_T = deepcopy(network.lines[:s_nom_opt])
        m = run_lopf(network, solver; rescaling=rescaling,formulation="angles_linear", investment_type="continuous")

        # write best results to network
        for l=1:nrow(network.lines)
            network.lines[:s_nom][l] = s_nom_0[l]
            network.lines[:s_nom_opt][l] = s_nom_opt_T[l]
            network.lines[:s_nom_extendable][l] = s_nom_extendable_0[l]
        end

    end

    # return model and iteration data
    return m, [objectives, capacities, reactances]

end