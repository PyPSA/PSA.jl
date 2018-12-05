using LightGraphs
using PyCall
using Missings
using DataFrames
using AxisArrays
const networkx = PyNULL()
copy!(networkx, pyimport("networkx" ))

include("axis_utils.jl"); include("graph.jl") 

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
    assert(in(typeof(snapshots), [Array{DateTime,1}, Array{Union{DateTime, Missings.Missing},1}]))

    snapshots = DateTime.(snapshots)
    n.snapshots = snapshots 
    for field=dynamic_components(n)
        for df_name=keys(getfield(n,field))
            if size(getfield(n,field)[df_name])[1]>0
                getfield(n,field)[df_name] = getfield(n,field)[df_name][snapshots, :]
            end
        end
    end
    if length(n.snapshot_weightings)>=0
        n.snapshot_weightings = n.snapshot_weightings[snapshots]
    end
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


function calculate_dependent_values!(n, f_o="1")
    function set_default(df, col, default)
        !in(col, df.axes[2].val) ? assign(df, fill(default, (size(df)[1])), col) : df
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

    # costs, minimize objective range by factor f_o 
    # f_o = float(f_o)
    # n.generators[:, "maintenance_cost"] .= n.generators[:, "maintenance_cost"]/f_o
    # n.generators[:, "marginal_cost"] .= n.generators[:, "marginal_cost"]/f_o
    # n.generators[:, "capital_cost"] .= n.generators[:, "capital_cost"]/f_o
    # n.links[:, "capital_cost"] .= n.links[:, "capital_cost"]/f_o
    # n.lines[:, "capital_cost"] .= n.lines[:, "capital_cost"]/f_o

    # check for infinity and replace for not having trouble in the lopf
    g_inf = n.generators[:, "p_nom_max"] .== Inf
    n.generators[g_inf, "p_nom_max"] = 1e9
    lin_inf = n.lines[:, "s_nom_max"] .== Inf
    n.lines[lin_inf, "s_nom_max"] = 1e9
    link_inf = n.links[:, "p_nom_max"] .== Inf
    n.links[link_inf, "p_nom_max"] = 1e9

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


"""
Get investment periods for infrastructure investments. Defaults to one
investment point at start of timeseries (normal lopf). If an investment period 
is given, one can decide to activate the first snapshot as an investment point 
or not. 
"""
function get_investment_periods(n, investment_period=nothing, 
                                invest_at_first_sn=false)

start = n.snapshots[1]
stop = n.snapshots[length(n.snapshots)]
periods = Dict("h" => Base.Dates.Hour(1), "d" => Base.Dates.Day(1),
"w" => Base.Dates.Week(1), "m" => Base.Dates.Month(1),
"y" => Base.Dates.Year(1))

    if isa(investment_period, Array)
        t = investment_period
    elseif investment_period in keys(periods)
        t = start:periods[investment_period]:stop
    else
        t = n.snapshots[1]:n.snapshots[1]
        investment_period == nothing ? nothing : warn("Not valid argument for investment period, falling back to first snapshot")
    end
    # array{DateTime,1} with snapshots when you invest + number of invesment times
    t = collect(t)
    if (length(t) > 1) & (invest_at_first_sn == false)
        t_ip = findin(n.snapshots, t)[2:end] 
        IP = length(t)-1    
    else
        t_ip = findin(n.snapshots, t) 
        IP = length(t) 
    end
    @assert length(t_ip) == IP 
    t_ip, IP
end


function aggregate_investments!(expr, start, var, t_ip)
    #set first snapshot
    if 1 ∈ t_ip
        expr[1,:] = start + var[1, :]
    else
        expr[1,:] = start'
    end

    #set following snapshots
    T = size(expr)[1]
    for t=2:T
        if t ∈ t_ip      
            ip = findin(t_ip, t)[1] 
            expr[t,:] = expr[t-1, :] + var[ip, :]
        else
            expr[t,:] = expr[t-1, :]
        end
    end
end

