using LightGraphs
using PyCall
using Missings
using DataFrames
using AxisArrays
const networkx = PyNULL()
copy!(networkx, pyimport("networkx" ))

"""
this function return all dynamic components of the network
"""
function dynamic_components(n)
    fields = String.(fieldnames(n))
    components = []; 
    for field=fields 
        field[end-1:end] == "_t" ? push!(components, field) : nothing  
    end
    Symbol.(components)
end

function static_components(n)
    fields = String.(fieldnames(n))
    deleteat!(fields, findin(fields, ["name"]))
    components = []
    for field=fields
        field[end-1:end] != "_t" ? push!(components, field) : nothing
    end
    Symbol.(components)
end


function set_snapshots!(n, snapshots)
    assert(typeof(snapshots) == Array{DateTime,1})
    first_sn, last_sn = snapshots[[1, end]]
    for field=dynamic_components(n)
        for df_name=keys(getfield(n,field))
            if size(getfield(n,field)[df_name])[1]>0
                getfield(n,field)[df_name] = getfield(n,field)[df_name][first_sn .. last_sn, :]
            end
        end
    end
    n.snapshots = snapshots
end


# align the column order of dynamic_component (e.g generators_t{["p_max_pu"]) with row order of 
# static component (generators[:name])
function align_component_order!(n)
    for comp = dynamic_components(n)
        order = getfield(n, Symbol(String(comp)[1:end-2])).axes[1]
        for attr in keys(getfield(n, comp))
            if ndims(getfield(n,comp)[attr]) > 1
                if size(getfield(n,comp)[attr])[2]== length(order) != 0
                    getfield(n,comp)[attr]= reindex(getfield(n,comp)[attr], columns = order)
                end
            end
        end
    end
end


function calculate_dependent_values!(n)
    function set_default(df, col, default)
        !in(col, df.axes[2].val) ? assign(df, fill(default, (size(df)[1]),1), col) : df
    end

    # snapshot_weighting
    (length(n.snapshot_weightings)==0?
        n.snapshot_weightings = AxisArray(fill(1, size(n.snapshots)[1]),Axis{:row}(n.snapshots)) : nothing)
    
        #buses
    defaults = [("v_nom", 1.)]
    for (col, default) in defaults
        n.buses = set_default(n.buses, col, default)
    end

    # generators
    defaults = [("p_nom_extendable", false), ("p_nom_max", Inf),("commitable", false),
                ("p_min_pu", 0), ("p_max_pu", 1), ("p_nom_min", 0),("capital_cost", 0),
                ("min_up_time", 0), ("min_down_time", 0), ("initial_status", true),
                ("p_nom", 0.),("marginal_cost", 0),("p_nom_opt", 0.), ("efficiency", 1.)]
    for (col, default) in defaults
        n.generators = set_default(n.generators, col, default)
    end

    # lines
    # map the bus v_nom to lines bus0
    vnom = []; for bus=string.(n.lines[:,"bus0"]) push!(vnom, n.buses[bus, "v_nom"]) end    
    n.lines = assign(n.lines, vnom, "v_nom")

    defaults = [("s_nom_extendable", false), ("s_nom_min", 0),("s_nom_max", Inf), ("s_nom", 0.),
                ("s_nom_min", 0), ("s_nom_max", Inf), ("capital_cost", 0), ("g", 0), ("s_nom_opt", 0.), 
                ("r", 0), ("b",0), ("g",0)]
    for (col, default) in defaults
        n.lines = set_default(n.lines, col, default)
    end
    n.lines = assign(n.lines, float.(n.lines[:,"x"])./float.((n.lines[:,"v_nom"]).^2), "x_pu")
    n.lines = assign(n.lines, float.(n.lines[:,"r"])./float.(n.lines[:,"v_nom"]).^2, "r_pu")
    n.lines = assign(n.lines, float.(n.lines[:,"b"]).*float.(n.lines[:,"v_nom"]).^2,   "b_pu")
    n.lines = assign(n.lines, float.(n.lines[:,"g"]).*float.(n.lines[:,"v_nom"]).^2,   "g_pu")

    # links
    defaults = [("p_nom_extendable", false), ("p_nom_max", Inf), ("p_min_pu", 0),
                ("p_max_pu", 1),("p_nom_min", 0), ("p_nom_max", Inf), ("capital_cost", 0),
                ("marginal_cost", 0), ("p_nom", 0.), ("efficiency", 1), ("p_nom_opt", 0.)]
    for (col, default) in defaults
        n.links = set_default(n.links, col, default)
    end

    # storage_units
    defaults = [("p_nom_min", 0), ("p_nom_max", Inf), ("p_min_pu", -1),
                ("p_max_pu", 1), ("marginal_cost", 0), ("efficiency_store", 1),
                ("cyclic_state_of_charge", false), ("p_nom_extendable", false),
                ("state_of_charge_initial", 0.), ("p_nom", 0.),("capital_cost", 0),
                ("efficiency_dispatch", 1), ("inflow", 0), ("p_nom_opt", 0.)]
    for (col, default) in defaults
        n.storage_units = set_default(n.storage_units, col, default)
    end

    # stores
    defaults = [("e_nom_min", 0), ("e_nom_max", Inf), ("e_min_pu", -1),
                    ("e_max_pu", 1), ("marginal_cost", 0), ("efficiency_store", 1),("capital_cost", 0),
                    ("efficiency_dispatch", 1),("inflow", 0), ("e_nom", 0.), ("e_nom_opt", 0.), 
                    ("max_hours", 0)]
    for (col, default) in defaults
        n.stores = set_default(n.stores, col, default)
    end

    # loads_t
    # could be directly done with assign
    for df_name=keys(n.loads_t)
        if size(n.loads_t[df_name])[1]>1
            for bus=n.loads.axes[1].val 
                if !in(bus, n.loads_t[df_name].axes[2].val)
                    n.loads_t[df_name] = set_default(n.loads_t[df_name], bus, 0)
                end
            end
        end
    end
    align_component_order!(n)
end


function get_switchable_as_dense(n, component, attribute, snapshots=0)
    """
    this function returns a matrix snapshot x component[attribute], so it collects time-varying and static attributes
    (e.g. "p_max_pu") of one component(e.g. "generators")
    """
    snapshots==0 ? snapshots = n.snapshots : nothing
    T, = size(snapshots)
    component_t = Symbol(component * "_t")
    component = Symbol(component)
    dense = getfield(n, component_t)[attribute] # time-varying data
    if size(dense)[1] > 0  # if there is a time-varying attribute
        static_value = getfield(n, component)[:,attribute]
        missing_cols = setdiff(static_value.axes[1].val, dense.axes[2].val)
        dense = assign(dense, repmat(float.(static_value[missing_cols].data)', T), missing_cols )
        reindex(dense, columns=static_value.axes[1])
    else
        static_value = getfield(n, component)[:,attribute]
        AxisArray(repmat(float.(static_value.data)', T), Axis{:snapshots}(n.snapshots), static_value.axes[1])
    end
end


function make_static_fallback_getter(df, static, selector=nothing)
    selector == nothing ? (selector = :) : nothing
    fullindex = static.axes[1].val[selector]
    static = static.data[selector]
    dynindex = df.axes[2].val
    (t,x)-> in(fullindex[x], dynindex) ? df[t, fullindex[x]] : static[x] 
    # end
end

function make_fallback_getter(df, def=1)
    (t,x)->in(x, names(df)) ? df[t, x] : def
end

function get_investment_periods(n, investment_period)
    # for investment_period
    start = n.snapshots[1]
    stop = n.snapshots[length(n.snapshots)]
    periods = Dict("h" => Base.Dates.Hour(1), "d" => Base.Dates.Day(1),
                   "w" => Base.Dates.Week(1), "m" => Base.Dates.Month(1),
                   "y" => Base.Dates.Year(1))

    if investment_period in keys(periods)
        t = start:periods[investment_period]:stop
    else
        t = n.snapshots[1]:n.snapshots[1]
        investment_period == nothing ? nothing : warn("Not valid argument for investment period, falling back to first snapshot")
    end
    # array{DateTime,1} with snapshots when you invest + number of invesment times
    collect(t)
    findin(n.snapshots, t), length(t) 
end


# -------------------------------------------------------------------------------------------------
# AxisArray functions

idx(dataframe) = Dict(zip(dataframe.axes[1].val, Iterators.countfrom(1)))
rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe.axes[1].val))

zsum(array) = length(array)>0 ? sum(array) : 0.
zdot(v1,v2) = length(v1)>0 ? dot(v1,v2) : 0.


# Use this funtion to assign a new (set of) column(s) or row(s) with given values.
function assign(df::AxisArray, values, index; axis=1)
    # This could also be done with AxisArrays.merge which breaks however with type 
    # Any.  
    in(typeof(index),[String, Symbol, DateTime]) ? index = [index] : nothing
    if axis == 1
        AxisArray([df.data values], 
                    axes(df)[1], 
                    Axis{axisnames(df)[2]}(append!(copy(axes(df)[2].val), index))) 
    elseif axis==2
        AxisArray([df.data; values'], 
                    Axis{axisnames(df)[1]}(append!(copy(axes(df)[1].val), index)),  
                    axes(df)[2])
    end
end

# might be a better way
# function assign(df::AxisArray, values, index; axis=2)
#     # This could also be done with AxisArrays.merge which breaks however with type 
#     # Any.  
#     in(typeof(index),[String, Symbol, DateTime]) ? index = [index] : nothing
#     axname = axisnames(df)
#     if axis == 1
#         values = AxisArray(values, Axis{axname[1]}(index), df.axes[2])
#     elseif axis==2
#         values = AxisArray(values, df.axes[1], Axis{axname[2]}(index))
#     end
#     cat(axis, df, values)
# end


# use this function to reorder, AxisArrays does not support this 
function reindex(df::AxisArray; index=nothing, columns=nothing)

    if index == nothing
        (newindex = :)
    else
        typeof(index) <: AxisArrays.Axis ? index = index.val : nothing 
        indexdict = Dict(zip(df.axes[1].val, Iterators.countfrom(1)))
        newindex = Int[]
        for i=index push!(newindex, getindex(indexdict, i)) end 
    end

    if columns == nothing
        (newcol = :)
    else
        typeof(columns) <: AxisArrays.Axis ? columns = columns.val : nothing 
        coldict = Dict(zip(df.axes[2].val, Iterators.countfrom(1)))
        newcol = Int[]
        for i=columns push!(newcol, getindex(coldict, i)) end 
    end

    df[newindex, newcol]
end

function grouped_array(A)
    groups = Dict()
    for element ∈ unique(A)
        groups[element] = findin(A, [element])
    end
    return groups
end

# like groupby in pandas 
function group(A, by_array, func; axis=1)
    by_array = grouped_array(by_array)
    grouped = []
    if ndims(A) == 1
        for index = keys(by_array)
            push!(grouped, func(A[by_array[index]]))
        end
    else
        if axis==1
            for index = keys(by_array)
                push!(grouped, func(A[by_array[index],:], axis))
            end
        elseif axis==2
            for index = keys(by_array)
                push!(grouped, func(A[:,by_array[index]], axis))
            end            
        end
        axname = axisnames(A)[axis]
        grouped = cat(axis, grouped...)
        axs = grouped.axes |> collect
        axs[axis] =  Axis{axname}(by_array |> keys |> collect .|> string)
        grouped = AxisArray(grouped.data, axs...) 
    end
    grouped
end


# --------------------------------------------------------------------------------------------------
# DataFrames functions


# auxilliaries
# idx(dataframe) = Dict(zip(dataframe[:name], Iterators.countfrom(1)))
# rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe[:name]))
idx_sym(dataframe) = Dict(zip(Symbol.(dataframe[:name]), Iterators.countfrom(1)))
rev_idx_sym(dataframe) = Dict(zip(Iterators.countfrom(1), Symbol.(dataframe[:name])))

# try to match pandas useful reindex with capavility of set_index, note that it also allows to index from 
# a duplicate axis. 
function reindex(df::DataFrame; index = nothing, columns = nothing, index_col=nothing, fill_value=missing)
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


function get_rows_of(df::DataFrame, names, index_col=:name)
    dict = Dict(zip(df[index_col], Iterators.countfrom(1)))
    positions = []
    for i=names push!(positions,getindex(idx_gen, i)) end 
    positions
end


function append_idx_col!(df::DataFrame)
    if typeof(df)==Vector{DataFrames.DataFrame}
        for df in df
            df[:idx] = collect(1:nrow(df))
        end
    else
        df[:idx] = collect(1:nrow(df))
    end
end

# -----------------------------------------------------------------------------------------------
# graph analysis

function to_graph(n)
    busidx = idx(n.buses)
    g = DiGraph(size(n.buses)[1])
    for l=1:size(n.lines)[1]
        add_edge!(g, busidx[n.lines[l,"bus0"]], busidx[n.lines[l,"bus1"]] )
    end
    for l=1:size(n.links)[1]
        add_edge!(g, busidx[n.links[l,"bus0"]], busidx[n.links[l,"bus1"]] )
    end
    return g
end

function incidence_matrix(n)
    busidx = idx(n.buses)
    lines = n.lines
    K = zeros(size(n.buses)[1],size(lines)[1])
    for l in 1:size(K)[2]
        K[busidx[lines[l,"bus0"]],l] = 1
        K[busidx[lines[l,"bus1"]],l] = -1
    end
    return K
end

function laplace_matrix(n)
    K = incidence_matrix(n)
    return K*K'
end

function ptdf_matrix(n)
    K = incidence_matrix(n)
    H = K' * pinv(K*K')
    return H .- H[:,1]
end

function get_cycles(n)
    busidx = idx(n.buses)
    g = networkx[:Graph]()
    g[:add_nodes_from](1:size(n.buses)[1])
    g[:add_edges_from]([(busidx[n.lines[l,"bus0"]], busidx[n.lines[l,"bus1"]]) for l=1:size(n.lines)[1]])
    networkx[:cycle_basis](g)
end

