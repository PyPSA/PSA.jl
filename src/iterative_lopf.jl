function iterative_lopf(network, solver, iterations; formulation::String="angles", objective::String="total", investment_type::String="continuous")

    # calculate_dependent_values!(network)
    x_0 = deepcopy(network.lines[:x])
    s_nom_0 = deepcopy(network.lines[:s_nom])

    objectives = Float64[]
    capacities = Array{Float64,1}[]
    reactances = Array{Float64,1}[]

    for k=1:iterations

        println("Run iteration $k")

        m = lopf(network, solver; formulation=formulation, objective=objective, investment_type=investment_type)

        push!(objectives, m.objVal)
        push!(capacities, network.lines[:s_nom_opt])
        push!(reactances, network.lines[:x])
        
        for l=1:nrow(network.lines)
            if network.lines[:s_nom_extendable][l]
                network.lines[:x][l] = ( s_nom_0[l]  * x_0[l] ) / network.lines[:s_nom_opt][l]
            end
        end
        
    end

    return [objectives, capacities, reactances]

end