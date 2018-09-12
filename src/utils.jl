using LightGraphs
using LightGraphs.SimpleGraphs

include("compat.jl")

function time_dependent_components(network)
    fields = String.(fieldnames(network))
    components = []
    for field=fields
        if field[end-1:end] == "_t"
            push!(components, field)
        end
    end
    Symbol.(components)
end


function static_components(network)
    fields = String.(fieldnames(network))
    components = []
    for field=fields
        if field[end-1:end] != "_t"
            push!(components, field)
        end
    end
    Symbol.(components)
end

function set_snapshots!(network, snapshots)
    for field=time_dependent_components(network)
        for df_name=keys(getfield(network,field))
            if nrow(getfield(network,field)[df_name])>0
                getfield(network,field)[df_name] = getfield(network,field)[df_name][snapshots, :]
            end
        end
    end
    network.snapshots = network.snapshots[snapshots,:]
end


function align_component_order!(network)
    components_t = time_dependent_components(network)
    for comp=components_t
        order = Symbol.(getfield(network, Symbol(String(comp)[1:end-2]))[:name])
        for attr in keys(getfield(network, comp))
            if length(getfield(network,comp)[attr])==length(order)
                getfield(network,comp)[attr]= getfield(network,comp)[attr][:, order]
            end
        end
    end
end

idx(dataframe) = Dict(zip(dataframe[:name], Iterators.countfrom(1)))
rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe[:name]))
idx_by(dataframe, col, values) = select_by(dataframe, col, values)[:idx]

function select_by(dataframe, col, selector)
    if length(findin(dataframe[col], selector))==0
        return dataframe[repeat(Bool[false],outer=nrow(dataframe)) , :]
    else
        mdict = Dict(zip(dataframe[col], Iterators.countfrom(1)))
        ids = Array{Int,1}(0)
        for i in selector
            push!(ids, mdict[i])
        end
        dataframe[ids,:]
    end
end
select_names(a, b) = select_by(a, :name, b)


function append_idx_col!(dataframe)
    if typeof(dataframe)==Vector{DataFrames.DataFrame}
        for df in dataframe
            df[:idx] = collect(1:nrow(df))
        end
    else
        dataframe[:idx] = collect(1:nrow(dataframe))
    end
end

function get_switchable_as_dense(network, component, attribute, snapshots=0)
    snapshots==0 ? snapshots = network.snapshots : nothing
    T = nrow(snapshots)
    component_t = Symbol(component * "_t")
    component = Symbol(component)
    dense = DataFrame()
    if in(attribute, keys(getfield(network, component_t)))
        dense = getfield(network, component_t)[attribute]
    end
    cols = Symbol.(getfield(network, component)[:name])
    not_included = String.(setdiff(cols, names(dense)))
    if length(not_included)>0
        attribute = Symbol.(attribute)
        df = select_names(getfield(network, component), not_included)
        df = names!(DataFrame(repmat(transpose(Array(df[attribute])), T)),
                Symbol.(not_included))
        dense = [dense df]
    end
    return dense[cols]
end

"""Returns a function which selects time dependent variables either
from the series network.component_t or from static valuenetwork.component"""
function select_time_dep(network, component, attribute; components=0)

    component_t = Symbol(component * "_t")
    component = Symbol(component)
    attribute_s = Symbol(attribute)

    df = getfield(network, component_t)[attribute]

    if components==0
        c_names = Symbol.(getfield(network, component)[:,:name])
        s = getfield(network, component)[:,attribute_s]
    else
        c_names = Symbol.(getfield(network, component)[components,:name])
        s = getfield(network, component)[components,attribute_s]
    end

    df_names = names(df)

    return (t,i) -> in(c_names[i], df_names) ? df[t, c_names[i]] : s[i]
end


function calculate_dependent_values!(network)
    function set_default(dataframe, col, default)
        !in(col, names(dataframe)) ? dataframe[col] = default : nothing
    end

    #buses
    defaults = [(:v_nom, 1.)]
    for (col, default) in defaults
        set_default(network.buses, col, default)
    end

    # generators
    defaults = [(:p_nom_extendable, false), (:p_nom_max, Inf),(:commitable, false),
                (:p_min_pu, 0), (:p_max_pu, 1), (:p_nom_min, 0),(:capital_cost, 0),
                (:min_up_time, 0), (:min_down_time, 0), (:initial_status, true),
                (:p_nom, 0.),(:marginal_cost, 0),(:p_nom_opt, 0.)]
    for (col, default) in defaults
        set_default(network.generators, col, default)
    end

    # lines
    network.lines[:v_nom]=select_names(network.buses, network.lines[:bus0])[:v_nom]
    defaults = [(:s_nom_extendable, true), (:s_nom_min, 0),(:s_nom_max, Inf), (:s_nom, 0.),
                (:capital_cost, 0), (:g, 0), (:s_nom_ext_min, 0.1), (:s_max_pu, 1.0),] 
    for (col, default) in defaults
        set_default(network.lines, col, default)
    end
    
    network.lines[:x_pu] = network.lines[:x]./(network.lines[:v_nom].^2)
    network.lines[:r_pu] = network.lines[:r]./(network.lines[:v_nom].^2)
    #network.lines[:b_pu] = network.lines[:b].*network.lines[:v_nom].^2
    #network.lines[:g_pu] = network.lines[:g].*network.lines[:v_nom].^2

    # links
    defaults = [(:p_nom_extendable, false), (:p_nom_max, Inf), (:p_min_pu, 0),
                (:p_max_pu, 1),(:p_nom_min, 0), (:p_nom_max, Inf), (:capital_cost, 0),
                (:marginal_cost, 0), (:p_nom, 0.), (:efficiency, 1)]
    for (col, default) in defaults
        set_default(network.links, col, default)
    end

    # storage_units
    defaults = [(:p_nom_min, 0), (:p_nom_max, Inf), (:p_min_pu, -1),
                (:p_max_pu, 1), (:marginal_cost, 0), (:efficiency_store, 1),
                (:cyclic_state_of_charge, false),
                (:state_of_charge_initial, 0.), (:p_nom, 0.),
                (:efficiency_dispatch, 1), (:inflow, 0)]
    for (col, default) in defaults
        set_default(network.storage_units, col, default)
    end

    # stores
    defaults = [(:e_nom_min, 0), (:e_nom_max, Inf), (:e_min_pu, -1),
                    (:e_max_pu, 1), (:marginal_cost, 0), (:efficiency_store, 1),
                    (:efficiency_dispatch, 1),(:inflow, 0), (:e_nom, 0.)]
    for (col, default) in defaults
        set_default(network.stores, col, default)
    end

    # loads_t
    for df_name=keys(network.loads_t)
        if nrow(network.loads_t[df_name])>1
            for bus=[bus for bus in network.loads[:name] if
                !in(Symbol(bus), names(network.loads_t[df_name]))]
                set_default(network.loads_t[df_name], bus, 0)
            end
        end
    end

end

function to_graph(network)
    busidx = idx(network.buses)
    g = DiGraph(length(busidx))
    for l in eachrow(network.lines)
        add_edge!(g, busidx[l[:bus0]], busidx[l[:bus1]] )
    end
    for l in eachrow(network.links)
        add_edge!(g, busidx[l[:bus0]], busidx[l[:bus1]] )
    end
    return g
end

function incidence_matrix(network)
    busidx = idx(network.buses)
    lines = network.lines
    K = zeros(nrow(network.buses),nrow(lines))
    for l in 1:size(K)[2]
        K[busidx[lines[l,:bus0]],l] = 1
        K[busidx[lines[l,:bus1]],l] = -1
    end
    return K
end

function laplace_matrix(network)
    K = incidence_matrix(network)
    return K*K'
end

function weighted_laplace_matrix(network)
    K = incidence_matrix(network)
    B = susceptance_matrix(network)
    return K*B*K'
end

function susceptance_matrix(network)
    return diagm(network.lines[:x].^(-1))
end

function ptdf_matrix(network)
    K = incidence_matrix(network)
    B = susceptance_matrix(network)
    L = weighted_laplace_matrix(network)
    H = B * K' * pinv(L)
    return H .- H[:,1]
end

function get_cycles(network)
    busidx = idx(network.buses)
    elist = [(busidx[l[:bus0]], busidx[l[:bus1]]) for l in eachrow(network.lines)]
    g = SimpleGraph(length(busidx))
    for e in elist add_edge!(g,e) end
    cycle_basis(g)
end

# following lift-and-project relaxation of Taylor2015a
function get_shortest_line_paths(network)
    busidx = idx(network.buses)

    # create elist only for lines with s_nom > 0
    elist = Tuple{Int64,Int64}[]
    for l in eachrow(network.lines)
        if l[:s_nom] > 0
            push!(elist, (busidx[l[:bus0]], busidx[l[:bus1]]))
        end
    end

    # create graph
    g = SimpleGraph(length(busidx))
    for e in elist add_edge!(g,e) end

    # set weights
    weights = zeros(nrow(network.lines), nrow(network.lines))
    for l in eachrow(network.lines)
        weights[busidx[l[:bus0]],busidx[l[:bus1]]] = l[:s_nom_step] * l[:x_step]
        weights[busidx[l[:bus1]],busidx[l[:bus0]]] = l[:s_nom_step] * l[:x_step]
    end

    line_paths = Array{Int64,1}[]
    for l in eachrow(network.lines)

        path_lg = a_star(g, busidx[l[:bus0]], busidx[l[:bus1]], weights)
        path = [(p.src,p.dst) for p in path_lg]
    
        lines_on_path = Int64[]
        for p in path
            for i=1:length(elist)
                if length(union(p, elist[i])) == length(p)
                    push!(lines_on_path, i)
                end
            end
        end

        push!(line_paths, lines_on_path)
    end

    return line_paths
end

function row_sum(df, row_id)
    if length(df[row_id,:]) == 0
        return 0.
    else
        return sum([df[row_id,i] for i in 1:length(df[row_id,:])])
    end
end

function get_iis(m::JuMP.Model)
    grb_model = MathProgBase.getrawsolver(internalmodel(m))
    num_constrs = Gurobi.num_constrs(grb_model)
    Gurobi.computeIIS(grb_model)
    iis_constrs = Gurobi.get_intattrarray(grb_model, "IISConstr",  1, num_constrs)
    m.linconstr[find(iis_constrs)]
end

# if slack is zero, relaxation would improve results
function get_slack(m::JuMP.Model)
    Gurobi.get_dblattrarray( m.internalModel.inner, "Slack", 1, Gurobi.num_constrs(m.internalModel.inner))
end

# active constraints are constraints with no slack
function get_active_constraints(m::JuMP.Model)

    slack = get_slack(m)
    
    function within_tolerance(x)
        (x <= 1e-9 && x >= -1e-9)
    end
    
    sort(collect(Set(find(within_tolerance,slack))))
end

function print_active_constraints!(m::JuMP.Model)

    slack = get_slack(m)
    active_constraints = get_active_constraints(m)

    for ac in active_constraints
        println(round(slack[ac],5),"\t\t",m.linconstr[ac])
    end

end

# TODO: choosing individual snapshots
function set_snapshots_start_end!(network, starting::String, ending::String)

    df = network.snapshots
    network.snapshots = df[(df[:name].<=ending).&(df[:name].>=starting),:]

    df = network.generators_t["p_max_pu"]
    network.generators_t["p_max_pu"] = df[(df[:name].<=ending).&(df[:name].>=starting),:]

    df = network.loads_t["p_set"]
    network.loads_t["p_set"] = df[(df[:name].<=ending).&(df[:name].>=starting),:]

end

# TODO: simple approach, more sophisticated method would cluster
function set_snapshots_sampling!(network, sampling_rate)

    T = nrow(network.snapshots)
    rows_to_delete = setdiff(collect(1:T),collect(1:sampling_rate:T))
    network.snapshots = deleterows!(network.snapshots,rows_to_delete)
    network.generators_t["p_max_pu"] = deleterows!(network.generators_t["p_max_pu"],rows_to_delete)
    network.loads_t["p_set"] = deleterows!(network.loads_t["p_set"],rows_to_delete)

end