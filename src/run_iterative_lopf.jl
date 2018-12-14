function run_iterative_lopf(network, solver, iterations; rescaling::Bool=false,formulation::String="angles_linear", investment_type::String="continuous", post_discretization::Bool=false, discretization_thresholds::Array{Float64,1}=[0.2,0.3], blockmodel=false, decomposition="benders")

    x_0 = deepcopy(network.lines[:x])
    s_nom_0 = deepcopy(network.lines[:s_nom])

    # debug
    println("x_0 = $x_0")
    println("s_nom_0 = $s_nom_0")

    objectives = Float64[]
    capacities = Array{Float64,1}[]
    reactances = Array{Float64,1}[]
    m = nothing

    for k=1:iterations

        println("Run iteration $k")

        m = run_lopf(network, solver; rescaling=rescaling,formulation=formulation, investment_type=investment_type, blockmodel=blockmodel, decomposition=decomposition)

        push!(objectives, m.objVal)
        push!(capacities, network.lines[:s_nom_opt])
        push!(reactances, network.lines[:x])

        # debug
        println("s_nom_opt_$k = $(network.lines[:s_nom_opt])")
        println("x_$k (before) = $(network.lines[:x])")

        for l=1:nrow(network.lines)
            if network.lines[:s_nom_extendable][l]
                if network.lines[:s_nom_opt][l] == 0.0
                    network.lines[:x][l] = 10e6 # reactance cannot take infinity values, instead choose prohibitively high value!
                else
                    network.lines[:x][l] = (x_0[l] * s_nom_0[l]) / network.lines[:s_nom_opt][l]
                end
            end
        end

        # debug
        println("x_$k (after) = $(network.lines[:x])")
            
    end

    # perform post discretization if selected
    if post_discretization

        # store the optimal solution of continuous optimisation
        s_nom_opt_continuous = deepcopy(network.lines[:s_nom_opt])
        s_nom_extendable_0 = deepcopy(network.lines[:s_nom_extendable])
        num_parallel_0 = deepcopy(network.lines[:num_parallel])
        m_opt = nothing
        threshold_opt = nothing
        threshold = discretization_thresholds[1]

        # iterate through all possible discretization thresholds
        println("LENGTH OF DISC THRESH IS $(length(discretization_thresholds))")
        if (length(discretization_thresholds) > 1)

            for threshold in discretization_thresholds

                println("#######")
                println("START EVALUATING THRESHOLD $threshold")
                println("#######")

                # round line extensions to integer
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

                # run lopf
                m_threshold = run_lopf(network, solver; rescaling=rescaling,formulation="angles_linear", investment_type="continuous")

                # compare to best solution in loop; better gets model
                if (threshold == discretization_thresholds[1]) || (m_threshold.objVal < m_opt.objVal)
                    println("#######")
                    println("CURRENT THRESHOLD BETTER OR FIRST ITERATION -- UPDATING threshold_opt to $threshold")
                    println("#######")
                    m_opt = m_threshold
                    threshold_opt = threshold
                end

                # reset network extendable circuits
                network.lines[:s_nom_extendable] = deepcopy(s_nom_extendable_0)

                println(network.lines[:s_nom_extendable])

            end
            
        end

        # run with optimal threshold choice
        println("#######")
        println("RUNNING AGAIN WITH OPTIMAL THHRESHOLD CHOICE $threshold_opt")
        println("#######")

        # round line extensions to integer
        for l=1:nrow(network.lines)

            if network.lines[:s_nom_extendable][l]

                num_parallel_extension = (s_nom_opt_continuous[l] / s_nom_0[l] - 1) * num_parallel_0[l]
                if mod(num_parallel_extension,1) >= threshold_opt
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