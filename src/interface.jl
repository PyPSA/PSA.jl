# imports
using PowerModels
pm = PowerModels

# replicates PowerModels network data for multiple snapshots without deepcopies
function replicate_light(sn_data::Dict{String,<:Any}, count::Int)
    @assert count > 0

    name = pm.get(sn_data, "name", "anonymous")

    mn_data = Dict{String,Any}(
        "nw" => Dict{String,Any}()
    )

    mn_data["multinetwork"] = true
    
    mn_data["name"] = "$(count) replicates of $(name)"

    for n in 1:count
        mn_data["nw"]["$n"] = sn_data
    end

    return mn_data
end

# Retrieves individual values of elements from PSA network
function val(element, attr::String, translator::Dict, busidx::Dict)

    # additional default values (mostly) not captured by PSA io
    defaults = Dict{String,Any}(
        "br_status" => 1,
        "tap" => 1.0,
        "shift" => 0.0,
        "transformer" => false,
        "b_fr" => 0.0,
        "b_to" => 0.0,
        "g_fr" => 0.0,
        "g_to" => 0.0,
        "angmin" => -pi/6,
        "angmax" => pi/6, 
        "vmin" => 0.9,
        "vmax" => 1.1
    )

    if in(translator[attr], keys(element.df.colindex))
        if occursin("bus",attr)
            return busidx[element[translator[attr]]]
        else
            return element[translator[attr]]
        end
    else
        if in(attr, keys(defaults))
            return defaults[attr]
        else
            return nothing
        end
    end
end

# adds elements of one component type to PowerModels network from PSA network
function get_component(network::PSA.network_mutable, comp::Symbol)

    # dictionary to translate variable names between PSA and PowerModels
    dictionary = Dict(
        :buses => Dict(
            "bus_i" => :name,
            "index" => :name,
            "base_kv" => :v_nom,
            "bus_type" => :control,
            "vmax" => :v_mag_pu_max,
            "vmin" => :v_mag_pu_min,
            "vm" => :v_mag_pu,
            "va" => :v_ang
            #"zone"
            #"area"
        ),
        :lines => Dict(
            "br_r" => :r_pu,
            "rate_a" => :s_nom,
            "shift" => nothing,# needs a default
            #"rate_b"
            "br_x" => :x_pu,
            #"rate_c"
            "g_fr" => :g_pu,
            "g_to" => :g_pu,
            "b_fr" => :b_pu,
            "b_to" => :b_pu,
            "f_bus" => :bus0,
            "br_status" => nothing, # needs a default
            "t_bus" => :bus1,
            "index" => :name,
            "angmin" => :v_ang_min,
            "angmax" => :v_ang_max,
            "transformer" => nothing, # needs a default
            "tap" => nothing # needs a default
        ),
    )

    attributes = dictionary[comp]
    idx = PSA.idx(getfield(network,comp))
    busidx = PSA.idx(getfield(network,:buses))
    component = Dict{String,Any}()
    for element=eachrow(getfield(network,comp))
        dict = Dict{String,Any}()
        for attr in keys(attributes)
            dict[attr] = val(element, attr, attributes, busidx)
        end
        component[repr(idx[element[:name]])] = dict
    end
    return component
end

# converts a PSA network to a PowerModels network (exclusively branch and bus data)
function convert_network(network::PSA.network_mutable; replicate::Bool=true)

    n = Dict{String, Any}()

    # read in relevant components from PSA
    components = [("bus",:buses), ("branch", :lines)]
    for (comp_pm, comp_psa) in components
       n[comp_pm] = get_component(network, comp_psa)
    end

    # must create dummy dictionaries for components which are not relevant
    dummy_components = ["dcline", "gen", "load", "shunt", "storage"]
    for comp in dummy_components
       n[comp] = Dict{String,Any}() 
    end

    # need to set reference bus
    n["bus"]["1"]["bus_type"] = 3 # if only one synchronous zone
    
    replicate ? n = replicate_light(n, nrow(network.snapshots)) : nothing
    
    return n
end

# single snapshot model recipe
function post_flow(gpm::GenericPowerModel)
    pm.variable_voltage(gpm)
    pm.variable_branch_flow(gpm, bounded=false)
    pm.constraint_voltage(gpm)
    for i in ids(gpm, :branch)
        pm.constraint_ohms_yt_from(gpm, i)
        pm.constraint_ohms_yt_to(gpm, i)
        pm.constraint_voltage_angle_difference(gpm, i)
    end
    for i in ids(gpm, :ref_buses)
        pm.constraint_theta_ref(gpm, i)
    end
end

# multiple snapshot model recipe
function post_mn_flow(gpm::GenericPowerModel)
    for (n, network) in nws(gpm)
        pm.variable_voltage(gpm, nw=n)
        pm.variable_branch_flow(gpm, nw=n, bounded=false)
        pm.constraint_voltage(gpm, nw=n)
        for i in ids(gpm, :branch, nw=n)
            pm.constraint_ohms_yt_from(gpm, i, nw=n)
            pm.constraint_ohms_yt_to(gpm, i, nw=n)
            pm.constraint_voltage_angle_difference(gpm, i, nw=n)
        end
        for i in ids(gpm, :ref_buses, nw=n)
            pm.constraint_theta_ref(gpm, i, nw=n)
        end
    end
end

# single snapshot model recipe for branch flow formulations
function post_flow_bf(gpm::GenericPowerModel)
    pm.variable_voltage(gpm)
    pm.variable_branch_flow(gpm, bounded=false)
    pm.variable_branch_current(gpm, bounded=false)

    for i in ids(gpm, :ref_buses)
        pm.constraint_theta_ref(gpm, i)
    end

    for i in ids(gpm, :branch)
        pm.constraint_flow_losses(gpm, i)
        pm.constraint_voltage_magnitude_difference(gpm, i)
        pm.constraint_branch_current(gpm, i)

        pm.constraint_voltage_angle_difference(gpm, i)
    end
end

# multiple snapshot model recipe for branch flow formulations
function post_mn_flow_bf(gpm::GenericPowerModel)
    for (n, network) in nws(gpm)
        pm.variable_voltage(gpm, nw=n)
        pm.variable_branch_flow(gpm, nw=n, bounded=false)
        pm.variable_branch_current(gpm, nw=n, bounded=false)

        for i in ids(gpm, :ref_buses, nw=n)
            pm.constraint_theta_ref(gpm, i, nw=n)
        end

        for i in ids(gpm, :branch, nw=n)
            pm.constraint_flow_losses(gpm, i, nw=n)
            pm.constraint_voltage_magnitude_difference(gpm, i, nw=n)
            pm.constraint_branch_current(gpm, i, nw=n)

            pm.constraint_voltage_angle_difference(gpm, i, nw=n)
        end
    end
end

""
function pm.variable_voltage(gpm::GenericPowerModel{T}; kwargs...) where T <: QCWRForm
    pm.variable_voltage_angle(gpm; kwargs...)
    pm.variable_voltage_magnitude(gpm; kwargs...)

    pm.variable_voltage_magnitude_sqr(gpm; kwargs...)
    pm.variable_voltage_product(gpm; kwargs...)

    pm.variable_voltage_angle_difference(gpm; kwargs...)
    pm.variable_voltage_magnitude_product(gpm; kwargs...)
    pm.variable_cosine(gpm; kwargs...)
    pm.variable_sine(gpm; kwargs...)
    variable_current_magnitude_sqr_unbounded(gpm; kwargs...)
end

""
function pm.variable_voltage(gpm::GenericPowerModel{T}; kwargs...) where T <: QCWRTriForm
    pm.variable_voltage_angle(gpm; kwargs...)
    pm.variable_voltage_magnitude(gpm; kwargs...)

    pm.variable_voltage_magnitude_sqr(gpm; kwargs...)
    pm.variable_voltage_product(gpm; kwargs...)

    pm.variable_voltage_angle_difference(gpm; kwargs...)
    pm.variable_voltage_magnitude_product(gpm; kwargs...)
    pm.variable_multipliers(gpm; kwargs...)
    pm.variable_cosine(gpm; kwargs...)
    pm.variable_sine(gpm; kwargs...)
    variable_current_magnitude_sqr_unbounded(gpm; kwargs...)
end

""
function variable_current_magnitude_sqr_unbounded(gpm::GenericPowerModel{T}; nw::Int=gpm.cnw, cnd::Int=gpm.ccnd) where T
    var(gpm, nw, cnd)[:cm] = @variable(gpm.model,
        [bp in ids(gpm, nw, :buspairs)], basename="$(nw)_$(cnd)_cm",
        lowerbound=0
    )
end

# converts PowerModel way of accessing variables to PSA way
function get_LN(network::PSA.network_mutable, powermodel::GenericPowerModel, variable::Symbol; reverse::Bool=false, ext::Bool=true)

    if typeof(powermodel) <: Union{DCPPowerModel, NFAPowerModel} && ( variable != :p || reverse )
        return Array{JuMP.Variable,2}
    end

    ln = network.lines
    T = nrow(network.snapshots)
    L = nrow(ln)
    busidx = PSA.idx(getfield(network,:buses))
    rev_lnidx = PSA.rev_idx(ln)
    ext_ln_b = ln[:s_nom_extendable]
    ext ? L_var = sum(ext_ln_b) : L_var = sum(.!ext_ln_b)
    vars = Array{Any,2}(L_var,T)
    reverse ? dir = (:bus1, :bus0) : dir = (:bus0, :bus1)
    l_var = 1
    for l=1:L
        if ext
            ln[ln[:name] .== rev_lnidx[l], :s_nom_extendable][1] == false ? continue : nothing
        else
            ln[ln[:name] .== rev_lnidx[l], :s_nom_extendable][1] == true ? continue : nothing
        end
        for t=1:T
            bus_fr = busidx[ln[ln[:name] .== rev_lnidx[l],dir[1]][1]]
            bus_to = busidx[ln[ln[:name] .== rev_lnidx[l],dir[2]][1]]
            vars[l_var,t] = var(powermodel, t, 1, variable)[(l,bus_fr,bus_to)]
        end
        l_var += 1
    end
    return convert(Array{JuMP.Variable,2}, vars)
end

""
function constraint_thermal_limit_from_ext(gpm::GenericPowerModel, i::Int, var::Union{JuMP.GenericAffExpr, JuMP.Variable}; nw::Int=gpm.cnw, cnd::Int=gpm.ccnd)
    if !haskey(con(gpm, nw, cnd), :sm_fr)
        con(gpm, nw, cnd)[:sm_fr] = Dict{Int,Any}() # note this can be a constraint or a variable bound
    end

    branch = ref(gpm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    pm.constraint_thermal_limit_from(gpm, nw, cnd, f_idx, var)
end


""
function constraint_thermal_limit_to_ext(gpm::GenericPowerModel, i::Int, var::Union{JuMP.GenericAffExpr, JuMP.Variable}; nw::Int=gpm.cnw, cnd::Int=gpm.ccnd)
    if !haskey(con(gpm, nw, cnd), :sm_to)
        con(gpm, nw, cnd)[:sm_to] = Dict{Int,Any}() # note this can be a constraint or a variable bound
    end

    branch = ref(gpm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    pm.constraint_thermal_limit_to(gpm, nw, cnd, t_idx, var)
end

""
function constraint_thermal_limit_from_fix(gpm::GenericPowerModel, i::Int, external_bound::Float64; nw::Int=gpm.cnw, cnd::Int=gpm.ccnd)
    if !haskey(con(gpm, nw, cnd), :sm_fr)
        con(gpm, nw, cnd)[:sm_fr] = Dict{Int,Any}() # note this can be a constraint or a variable bound
    end

    branch = ref(gpm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)

    pm.constraint_thermal_limit_from(gpm, nw, cnd, f_idx, external_bound)
end


""
function constraint_thermal_limit_to_fix(gpm::GenericPowerModel, i::Int, external_bound::Float64; nw::Int=gpm.cnw, cnd::Int=gpm.ccnd)
    if !haskey(con(gpm, nw, cnd), :sm_to)
        con(gpm, nw, cnd)[:sm_to] = Dict{Int,Any}() # note this can be a constraint or a variable bound
    end

    branch = ref(gpm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    t_idx = (i, t_bus, f_bus)

    pm.constraint_thermal_limit_to(gpm, nw, cnd, t_idx, external_bound)
end