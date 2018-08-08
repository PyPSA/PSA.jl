function run_iterative_lopf(network, solver, iterations; formulation::String="angles", objective::String="total", investment_type::String="continuous")

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

        m = run_lopf(network, solver; formulation=formulation, objective=objective, investment_type=investment_type)

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

    return m, [objectives, capacities, reactances]

end