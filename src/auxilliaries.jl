

function time_dependent_components(network)
    return [field for field=fieldnames(network) if String(field)[end-1:end]=="_t"]
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

# auxilliary functions
idx(dataframe) = Dict(zip(dataframe[:name], Iterators.countfrom(1)))
rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe[:name]))
select_names(dataframe, names) = dataframe[findin(dataframe[:name], names),:]
select_by(dataframe, col, values) = dataframe[findin(dataframe[col], values),:]
idx_by(dataframe, col, values) = select_by(dataframe, col, values)[:idx]

function to_symbol(str)
    if typeof(str)==String
        return Symbol(replace(str, " ", "_"))
    elseif (typeof(str)==Vector{String} || typeof(str) == DataArrays.DataArray{String,1})
        return [Symbol(replace(x, " ", "_")) for x in str]
    elseif typeof(str)==Int
        return Symbol("$str")
    else
        return str
    end
end

function to_string(sym)
    if typeof(sym)==Symbol
        return replace(String(sym), "_", " ")
    elseif (typeof(sym)==Vector{Symbol} || typeof(sym) == DataArrays.DataArray{String,1})
        return [replace(String(x), "_", " ") for x in sym]
    end
end

function append_idx_col!(dataframe)
    if typeof(dataframe)==Vector{DataFrames.DataFrame}
        for df in dataframe
            df[:idx] = 1:nrow(df)
        end
    else
        dataframe[:idx] = 1:nrow(dataframe)
    end
end

function get_switchable_as_dense(network, component, attribute, snapshots=0)
    snapshots==0 ? snapshots = network.snapshots : nothing
    T = nrow(snapshots)
    component_t = to_symbol(component * "_t")
    component = to_symbol(component)
    dense = DataFrame()
    if in(attribute, keys(getfield(network, component_t)))
        dense = getfield(network, component_t)[attribute]
    end
    cols = to_symbol(collect(getfield(network, component)[:name]))
    not_included = [to_string(c) for c=cols if !in(c,names(dense))]
    if length(not_included)>0
        attribute = to_symbol(attribute)
        df = select_names(getfield(network, component), not_included)
        df = names!(DataFrame(repmat(transpose(Array(df[attribute])), T)),
                to_symbol(not_included))
        dense = [dense df]
    end
    return dense[cols]
end


function calculate_dependent_values!(network)
    function set_default(dataframe, col, default)
        !in(col, names(dataframe)) ? dataframe[col] = default : nothing
    end
    defaults = [(:p_nom_extendable, false), (:p_nom_max, Inf),(:commitable, false),
                (:p_min_pu, 0), (:p_max_pu, 1), (:p_nom_min, 0),(:capital_cost, 0),
                (:min_up_time, 0), (:min_down_time, 0), (:initial_status, true),
                (:p_nom, NaN)]
    for (col, default) in defaults
        set_default(network.generators, col, default)
    end
    defaults = [(:s_nom_extendable, false), (:s_nom_min, 0),(:s_nom_max, Inf),
                (:s_nom_min, 0), (:s_nom_max, Inf), (:capital_cost, 0)]
    for (col, default) in defaults
        set_default(network.lines, col, default)
    end
    defaults = [(:p_nom_extendable, false), (:p_nom_max, Inf), (:p_min_pu, 0),
                (:p_max_pu, 1),(:p_nom_min, 0), (:p_nom_max, Inf), (:capital_cost, 0)]
    for (col, default) in defaults
        set_default(network.links, col, default)
    end
        defaults = [(:p_nom_min, 0), (:p_nom_max, Inf), (:p_min_pu, -1),
                    (:p_max_pu, 1), (:marginal_cost, 0), (:efficiency_store, 1),
                    (:efficiency_dispatch, 1)]
        for (col, default) in defaults
            set_default(network.storage_units, col, default)
    end
    for df_name=keys(network.loads_t)
        if nrow(network.loads_t[df_name])>1
            for bus=[bus for bus in network.buses[:name] if 
                !in(to_symbol(bus), names(network.loads_t[df_name]))]
                
            end
        end
    end 
end
