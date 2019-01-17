using Missings, DataFrames, AxisArrays, NamedTuples
const axes = Base.axes

include("axis_utils.jl"); include("graph.jl") 

function components(n)
    n |> typeof |> fieldnames .|> String |> collect
end

"""
Returns all dynamic components of the network
"""
function dynamic_components(n)
    dyns = []; 
    for c = components(n) 
        c[end-1:end] == "_t" ? push!(dyns, c) : nothing  
    end
    Symbol.(dyns)
end

function static_components(n)
    comps = components(n)
    filter!(c -> c != "name", comps)
    stats = []
    for c = comps
        c[end-1:end] != "_t" ? push!(stats, c) : nothing
    end
    Symbol.(stats)
end


"""
Replacing function, replaces attribute attr of a time-dependant component comp 
in network n by the value val. 
"""
function replace_attribute!(n, comp, attr, val)
    comp, attr = Symbol.([comp, attr])
    setindex(getfield(n, comp), attr, val) |> nt -> setfield!(n, comp, nt)
end

function set_snapshots!(n, snapshots)
    @assert(in(typeof(snapshots), [Array{DateTime,1}, 
            Array{Union{DateTime, Missings.Missing},1}]))

    snapshots = DateTime.(snapshots)
    n.snapshots = snapshots 
    for c = dynamic_components(n)
        for attr = keys(getfield(n, c))
            if size(getfield(n, c)[attr])[1] > 0
                val = getfield(n,c)[attr][snapshots, :]
                replace_attribute!(n, c, attr, val)
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
    for c = dynamic_components(n)
        order = getfield(n, Symbol(String(c)[1:end-2])).axes[1]
        for attr in keys(getfield(n, c))
            if ndims(getfield(n,c)[attr]) > 1
                if size(getfield(n,c)[attr])[2]== length(order) != 0
                    val = reindex(getfield(n,c)[attr], columns = order)
                    replace_attribute!(n, c, attr, val)
                end
            end
        end
    end
end


function calculate_dependent_values!(n)
    function set_default(df, col, default)
        !in(col, df.axes[2].val) ? assign(df, fill(default, (size(df)[1])), col) : df
    end

    # snapshot_weighting
    (length(n.snapshot_weightings)==0 ?
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
    "v_nom" ∈ n.lines.axes[2].val && sum(ismissing.(n.lines[:, "v_nom"])) == 0 ? nothing : (
        vnom = []; for bus=string.(n.lines[:,"bus0"]) push!(vnom, n.buses[bus, "v_nom"]) end;    
        n.lines = assign(n.lines, vnom, "v_nom") )

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
                ("efficiency_dispatch", 1), ("inflow", 0), ("p_nom_opt", 0.), ("max_hours", 1) ]
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

    # check for infinity and replace for not having trouble in the lopf
    for comp in static_components(n)
        data = getfield(n, comp)
        if sum(ismissing.(data)) != 0 || sum(skipmissing(BitArray(data.==NaN))) != 0
            info("Component $comp has $(sum(ismissing.(data))) missing values")
        end
        size(data)[1] > 0? data[(.!ismissing.(data)) .& (data .== Inf)] = 1e7 : nothing
    end

end

function scale_cost!(n, f_o="1")
"""""
this function scales the costs by the factor f_0 to make the matrix range smaller in the lopf
    
"""""
    cost = ["maintenance_cost", "capital_cost", "marginal_cost"]
        for comp in static_components(n)
            data = getfield(n, comp)
            if size(data)[1] > 0 && comp != :snapshots && comp != :snapshot_weightings
                for c in cost
                    c in data.axes[2][:]? data[:, c] .= data[:,c]/float(f_o) : nothing
                end
            end
        end
end



"""
this function returns a matrix snapshot x component[attribute], so it collects time-varying and static attributes
(e.g. "p_max_pu") of one component(e.g. "generators")
"""
function get_switchable_as_dense(n, component, attribute, snapshots=0)
    snapshots==0 ? snapshots = n.snapshots : nothing
    T, = size(snapshots)
    c_t = Symbol(component * "_t")
    c = Symbol(component)
    dense = getfield(n, c_t)[Symbol(attribute)] # dynamic data
    if size(dense)[1] > 0  # if dynamic data not empty
        static = getfield(n, c)[:, attribute]
        missing_cols = setdiff(static.axes[1].val, dense.axes[2].val)
        dense = assign(dense, repeat(float.(static[missing_cols].data)', T), missing_cols)
        reindex(dense, columns=static.axes[1])
    else
        static = getfield(n, c)[:, attribute]
        AxisArray(repeat(float.(static.data)', T), 
            Axis{:snapshots}(n.snapshots), static.axes[1])
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
periods = Dict("h" => Hour(1), "d" => Day(1), "w" => Week(1), 
               "m" => Month(1), "y" => Year(1))

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


function aggregate_investments!(expr, start, var, t_ip, invest_at_first_sn)
    #set first snapshot
    if 1 ∈ t_ip # || invest_at_first_sn==true
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
    @show(count(iszero, expr))
end
 # ------------------------------------------------------------------
function get_unused_capacity(n, dict, comp, attribute1, attribute2)
    (data_opt, data_flow) = (getfield(n, comp)[attribute1], getfield(n,comp)[attribute2])
    data_flow = data_flow[data_opt.axes[2]]   # get just the extendables
    num = count(sum(round.(data_opt),1) .== 0)
    count(sum(round.(data_opt),1) .== 0) > 0? info("in $comp there are $num unused") : nothing
    push!(dict, string(comp, "_unused") => (data_opt - abs.(data_flow)))
end

function get_summary(n)
    summary = Dict()
    #get unused capacities
    get_unused_capacity(n, summary, :generators_t, "p_nom_opt", "p")
    get_unused_capacity(n, summary, :lines_t, "s_nom_opt", "p0")
    get_unused_capacity(n, summary, :links_t, "p_nom_opt", "p0")
    
    # get p_nom_opt after carrier
    dict = Dict(n.generators.axes[1][i] => n.generators[i, "super_carrier"] 
            for i = 1:(size(n.generators)[1]))
    axes = unique(collect(n.generators[:, "super_carrier"]))
    

    attribute = ["p_nom_opt"]
    for comp in attribute
        list_carrier = fill(0.0, (length(n.snapshots), length(axes)))
        sum = (AxisArray(list_carrier, Axis{:time}(n.snapshots), 
                                    Axis{:col}(axes)))
        for sn=1:length(n.snapshots)
            t = n.generators_t[comp].axes[1][sn]
            for g=1:(size(n.generators_t[comp])[2])
                generator = n.generators_t[comp].axes[2][g]
                value = n.generators_t[comp][t,g]
                carrier = dict[generator]
                sum[t, carrier] += value 
            end
        end
        push!(summary, string(comp, "_carrier") => sum)
    end
    
    # check if emergeny generators are needed
    em_gen= n.generators[:, "super_carrier"] .== "emergency"
    em_p_nom_opt = round.(n.generators_t["p_nom_opt"][1, em_gen])
    indx = find(em_p_nom_opt .!= 0)
    length(indx)>0? @show(n.generators[em_gen, :][indx, :], n.generators_t["p_nom_opt"][:, em_gen][:,indx]) : nothing

    return summary
end

