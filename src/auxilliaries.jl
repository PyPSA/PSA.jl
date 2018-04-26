using LightGraphs
using PyCall
using Missings
using DataFrames
const networkx = PyNULL()
copy!(networkx, pyimport("networkx" ))


function dynamic_components(network)
    fields = String.(fieldnames(network))
    components = []; 
    for field=fields 
        field[end-1:end] == "_t" ? push!(components, field) : nothing  
    end
    Symbol.(components)
end

function static_components(network)
    fields = String.(fieldnames(network))
    components = []
    for field=fields
        field[end-1:end] != "_t" ? push!(components, field) : nothing
    end
    Symbol.(components)
end

function set_snapshots!(network, snapshots)
    for field=dynamic_components(network)
        for df_name=keys(getfield(network,field))
            if nrow(getfield(network,field)[df_name])>0
                getfield(network,field)[df_name] = getfield(network,field)[df_name][snapshots, :]
            end
        end
    end
    network.snapshots = network.snapshots[snapshots,:]
end


function align_component_order!(network)
    # align the column order of dynamic_component (e.g generators_t{["p_max_pu"]) with row order of 
    # static component (generators[:name])
    for comp = dynamic_components(network)
        order = Symbol.(getfield(network, Symbol(String(comp)[1:end-2]))[:name])
        for attr in keys(getfield(network, comp))
            if length(getfield(network,comp)[attr])==length(order)
                getfield(network,comp)[attr]= getfield(network,comp)[attr][:, order]
            end
        end
    end
end

# auxilliary funcitons
idx(dataframe) = Dict(zip(dataframe[:name], Iterators.countfrom(1)))
rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe[:name]))


# try to match pandas useful reindex with capavility of set_index, note that it also allows to index from 
# a duplicate axis. 
function reindex(df; index = nothing, columns = nothing, index_col=nothing, fill_value=missing)
    df = copy(df)
    if columns != nothing
        columns = Symbol.(columns)        
        for col=setdiff(columns, names(df))  df[col] = fill_value  end
        ordercols = Int64[]; for col=columns push!(ordercols, df.colindex[col]) end
        columns = names(df)[ordercols]
    else
        columns=names(df)
    end

    # sort columns such that index_col if existent or needed is at first place
    if index != nothing 
        if index_col == nothing
            index_col=:name
        end
        in(index_col, columns) ? deleteat!(columns, findin(columns, [index_col])) : nothing
        columns = append!([index_col], columns)        
    end
    df = df[columns]
    
    if index != nothing
        # deal with nonincluded index values
        not_in_index_col = setdiff(Array(index), df[index_col] )
        if length(not_in_index_col) > 0
            missing_rows = DataFrame( [not_in_index_col], [index_col])
            for col=columns[2:end] missing_rows[col] = missing end
            df = vcat(df, missing_rows)
        end
        
        # reorder aligned to index, first check if index_col is unique, if not takes go through every entry
        if any(nonunique(df[[index_col]]))
            warn("Indexing from a duplicated axis")
            orderindex = Int64[]; for i=index append!(orderindex , findin((df[index_col]), [i]) ) end
        else
            dict = Dict(zip(df[index_col], Iterators.countfrom(1)))
            orderindex = Int64[]; for i=index push!(orderindex, dict[i]) end
        end
       df[orderindex, :]
    else 
        df
    end
end



# function select_by(dataframe, col, selector)
#     if length(findin(dataframe[col], selector))==0
#         return dataframe[repeat(Bool[false],outer=nrow(dataframe)) , :]
#     else
#         mdict = Dict(zip(dataframe[col], Iterators.countfrom(1)))
#         ids = Array{Int,1}(0)
#         for i in selector
#             push!(ids, mdict[i])
#         end
#         dataframe[ids,:]
#     end
# end
# select_names(a, b) = select_by(a, :name, b)


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
        df = reindex(getfield(network, component), index=not_included)
        df = names!(DataFrame(repmat(transpose(Array(df[attribute])), T)),
                Symbol.(not_included))
        dense = [dense df]
    end
    return dense[cols]
end


function make_fallback_getter(df)
    def = 1.
    (t,x)->in(x, names(df)) ? df[t, x] : def
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
    network.lines[:v_nom]=reindex(network.buses, index=network.lines[:bus0])[:v_nom]
    defaults = [(:s_nom_extendable, false), (:s_nom_min, 0),(:s_nom_max, Inf), (:s_nom, 0.),
                (:s_nom_min, 0), (:s_nom_max, Inf), (:capital_cost, 0), (:g, 0)]
    for (col, default) in defaults
        set_default(network.lines, col, default)
    end
    network.lines[:x_pu] = network.lines[:x]./(network.lines[:v_nom].^2)
    network.lines[:r_pu] = network.lines[:r]./(network.lines[:v_nom].^2)
    network.lines[:b_pu] = network.lines[:b].*network.lines[:v_nom].^2
    network.lines[:g_pu] = network.lines[:g].*network.lines[:v_nom].^2

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

# function to_graphx(network)
#     busidx = idx(network.buses)
#     g = networkx[:Graph]()
#     g[:add_nodes_from](busidx)
#     g[:add_edges_from]([(busidx[l[:bus0]], busidx[l[:bus1]]) for l in eachrow(network.lines)])
#     return g
# end

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

function ptdf_matrix(network)
    K = incidence_matrix(network)
    H = K' * pinv(K*K')
    return H .- H[:,1]
end

function get_cycles(network)
    busidx = idx(network.buses)
    g = networkx[:Graph]()
    g[:add_nodes_from](busidx)
    g[:add_edges_from]([(busidx[l[:bus0]], busidx[l[:bus1]]) for l in eachrow(network.lines)])
    networkx[:cycle_basis](g)
end


function row_sum(df, row_id)
    if length(df[row_id,:]) == 0
        return 0.
    else
        return sum([df[row_id,i] for i in 1:length(df[row_id,:])])
    end
end
