
using AxisArrays, NCDatasets, NetCDF
include("/home/vres/data/lisa/jl/PSA.jl/src/auxilliaries.jl")

function char_to_string(data)
    # data = copy(data)
    data = Array{UInt8}(data)
    squeeze(mapslices(i -> unsafe_string(pointer(i)), data, 1), 1)
end



function reformat(data)
    if typeof(data) == Array{Int8,1}
        Array{Bool}(data)
    elseif typeof(data) == Array{UInt8,2}
        char_to_string(data)
    elseif typeof(data) == Array{Int64, 1}
        Array{Float64}(data)
    else
        data    
    end
end




# function nc_format(data::NCDatasets.CFVariable)
#     display(data.var)
#     if typeof(data) == NCDatasets.CFVariable{String,String,1}
#         Array{String}(data) 
#     elseif typeof(data) == NCDatasets.CFVariable{Char,Char,2}
#         replace.(prod(string.(data), 1), "\0", "")[1,:]
#     elseif typeof(data) == NCDatasets.CFVariable{Int64,Float64,1}
#         float(data)
#     elseif typeof(data) == NCDatasets.CFVariable{Int8,Float64,1}
#     # elseif in(:attrib, fieldnames(data))
#     #     in(data.attrib, "dtype") & (data.attrib["dtype"] == "bool") ? Array{Bool}(data) : data
#         Array{Bool}(data)
#     else
#         data
#     end
# end


function import_from_netcdf(path)
    n = axNetwork()
    ds_keys = keys(NCDatasets.Dataset(path)) 
    components = static_components(n)
    for comp = string.(components)
        display(comp)
        if any(contains.(ds_keys, comp))
            if comp == "snapshots"
                timeattr = split(ncgetatt(path, comp, "units"), " ")
                sns = DateTime(timeattr[3])+Dates.Hour.(ncread(path, comp))
                weightings = ncread(path, "snapshots_weightings")
                @show(weightings)
                setfield!(n, Symbol(comp),
                        AxisArray(  [sns weightings],
                        Axis{:index}(1:length(sns)), Axis{:properties}([comp, "weightings"]) ))             
            else
                index = reformat(ncread(path, comp * "_i"))
                props = ds_keys[find(contains.(ds_keys,   Regex("$(comp)_(?!(i\$|t_))"))) ]
                cols = props
                setfield!(n, Symbol(comp),
                        AxisArray(  hcat([reformat(ncread(path, key)) for key = props]...) ,
                        Axis{:index}(index), Axis{:properties}(replace.(props, comp * "_", ""))) )
            end
        end
    end
    components_t = dynamic_components(n)
    for comp=components_t
        for attr in keys(getfield(n, comp))
            comp_stat = Symbol(String(comp)[1:end-2])
        
            if in("$(comp)_$attr", ds_keys)
                getfield(n,comp)[attr]= (
                        AxisArray( ncread(path , "$(comp)_$attr")',
                                Axis{:snapshots}(sns), 
                                Axis{comp_stat}(reformat(ncread(path, "$(comp)_$(attr)_i")) ) )             
                )
            end
        end
    end
    n;
end




mutable struct network_mutable
    buses::AxisArray
    generators::AxisArray
    loads::AxisArray
    lines::AxisArray
    links::AxisArray
    storage_units::AxisArray
    stores::AxisArray
    transformers::AxisArray
    carriers::AxisArray
    global_constraints::AxisArray
    buses_t::Dict{String,AxisArray}
    generators_t::Dict{String,AxisArray}
    loads_t::Dict{String,AxisArray}
    lines_t::Dict{String,AxisArray}
    links_t::Dict{String,AxisArray}
    storage_units_t::Dict{String,AxisArray}
    stores_t::Dict{String,AxisArray}
    transformers_t::Dict{String,AxisArray}
    snapshots::AxisArray
end

function axNetwork(
# static
    buses = AxisArray(repeat([Any[]], outer=16),
            [:name, :v_nom, :type, :x, :y, :carrier, :v_mag_pu_set, :v_mag_pu_min,
            :v_mag_pu_max, :control, :sub_network, :p, :q,
            :v_mag_pu, :v_ang, :marginal_price]),

    carriers = AxisArray(repeat([Any[]], outer=2),
                 [:name, :co2_emissions]),

    generators = AxisArray(repeat([Any[]], outer=31),
            [:name, :bus, :control, :type, :p_nom, :p_nom_extendable,
            :p_nom_min, :p_nom_max, :p_min_pu, :p_max_pu, :p_set, :q_set,
            :sign, :carrier, :marginal_cost, :capital_cost, :efficiency,
            :committable, :start_up_cost, :shut_down_cost, :min_up_time,
            :min_down_time, :initial_status, :ramp_limit_up, :ramp_limit_down,
            :ramp_limit_start_up, :ramp_limit_shut_down, :p, :q, :p_nom_opt, :status]),

    global_constraints = AxisArray(repeat([Any[]], outer=6),
            [:name, :type, :carrier_attribute, :sense, :constant, :mu]),

    line_types = AxisArray(repeat([Any[]], outer=9),
            [:name, :f_nom, :r_per_length, :x_per_length, :c_per_length,
            :i_nom, :mounting, :cross_section, :references]),

    lines = AxisArray(repeat([Any[]], outer=33),
            [:name, :bus0, :bus1, :type, :x, :r, :g, :b,
            :s_nom, :s_nom_extendable, :s_nom_min, :s_nom_max,
            :s_max_pu, :capital_cost, :length, :terrain_factor, :num_parallel,
            :v_ang_min, :v_ang_max, :sub_network, :p0, :q0, :p1, :q1,
            :x_pu, :r_pu, :g_pu, :b_pu, :x_pu_eff, :r_pu_eff,
            :s_nom_opt, :mu_lower, :mu_upper]),

    links = AxisArray(repeat([Any[]], outer=21),
            [:name, :bus0, :bus1, :type, :efficiency, :p_nom, :p_nom_extendable,
            :p_nom_min, :p_nom_max, :p_set, :p_min_pu, :p_max_pu, :capital_cost,
            :marginal_cost, :length, :terrain_factor, :p0, :p1,
            :p_nom_opt, :mu_lower, :mu_upper]),

    loads = AxisArray(repeat([Any[]], outer=8),
            [:name, :bus, :type, :p_set, :q_set, :sign, :p, :q]),

    shunt_impendances = AxisArray(repeat([Any[]], outer=9),
            [:name, :bus, :g, :b, :sign, :p, :q, :g_pu, :b_pu]),

    storage_units = AxisArray(repeat([Any[]], outer=29),
            [:name, :bus, :control, :type, :p_nom, :p_nom_extendable, :p_nom_min,
            :p_nom_max, :p_min_pu, :p_max_pu, :p_set, :q_set, :sign, :carrier,
            :marginal_cost, :capital_cost, :state_of_charge_initial, :state_of_charge_set,
            :cyclic_state_of_charge, :max_hours, :efficiency_store, :efficiency_dispatch,
            :standing_loss, :inflow, :p, :q, :state_of_charge, :spill, :p_nom_opt]),

    stores = AxisArray(repeat([Any[]], outer=22),
            [:name, :bus, :type, :e_nom, :e_nom_extendable, :e_nom_min,
            :e_nom_max, :e_min_pu, :e_max_pu, :e_initial, :e_cyclic, :p_set,
            :cyclic_state_of_charge, :q_set, :sign, :marginal_cost,
            :capital_cost, :standing_loss, :p, :q, :e, :e_nom_opt]),

    transformer_types = AxisArray(repeat([Any[]], outer=16),
            [:name, :f_nom, :s_nom, :v_nom_0, :v_nom_1, :vsc, :vscr, :pfe,
            :i0, :phase_shift, :tap_side, :tap_neutral, :tap_min, :tap_max, :tap_step, :references]),

    transformers = AxisArray(repeat([Any[]], outer=36),
            [:name, :bus0, :bus1, :type, :model, :x, :r, :g, :b, :s_nom,
            :s_nom_extendable, :s_nom_min, :s_nom_max, :s_max_pu, :capital_cost,
            :num_parallel, :tap_ratio, :tap_side, :tap_position, :phase_shift,
            :v_ang_min, :v_ang_max, :sub_network, :p0, :q0, :p1, :q1, :x_pu,
            :r_pu, :g_pu, :b_pu, :x_pu_eff, :r_pu_eff, :s_nom_opt, :mu_lower, :mu_upper]),

# time_dependent
    buses_t=Dict{String,AxisArray}([("marginal_price",AxisArray([])), ("v_ang", AxisArray([])),
            ("v_mag_pu_set",AxisArray([])), ("q", AxisArray([])),
            ("v_mag_pu", AxisArray([])), ("p", AxisArray([])), ("p_max_pu", AxisArray([]))]),
    generators_t=Dict{String,AxisArray}(
            [("marginal_price",AxisArray([])), ("v_ang", AxisArray([])),
            ("v_mag_pu_set",AxisArray([])), ("q", AxisArray([])),
            ("v_mag_pu", AxisArray([])), ("p", AxisArray([])),
            ("p_min_pu", AxisArray([])), ("p_max_pu", AxisArray([]))
            ]),
    loads_t=Dict{String,AxisArray}([
            ("q_set",AxisArray([])), ("p_set", AxisArray([])),
            ("q", AxisArray([])), ("p", AxisArray([]))
            ]),
    lines_t=Dict{String,AxisArray}([("q0",AxisArray([])), ("q1", AxisArray([])),
            ("p0",AxisArray([])), ("p1", AxisArray([])),
            ("mu_lower", AxisArray([])), ("mu_upper", AxisArray([])),
            ("s_nom_opt", AxisArray([]))]),
    links_t=Dict{String,AxisArray}([("p_min_pu",AxisArray([])), ("p_max_pu", AxisArray([])),
            ("p0",AxisArray([])), ("p1", AxisArray([])),
            ("mu_lower", AxisArray([])), ("mu_upper", AxisArray([])),
            ("efficiency",AxisArray([])), ("p_set", AxisArray([])),
            ("p_nom_opt", AxisArray([]))]),

    storage_units_t=Dict{String,AxisArray}([("p_min_pu",AxisArray([])), ("p_max_pu", AxisArray([])),
        ("inflow",AxisArray([])), ("mu_lower", AxisArray([])), ("mu_upper", AxisArray([])),
        ("efficiency",AxisArray([])), ("p_set", AxisArray([]))]),
    stores_t=Dict{String,AxisArray}([("e_min_pu",AxisArray([])), ("e_max_pu", AxisArray([])),
        ("inflow",AxisArray([])), ("mu_lower", AxisArray([])), ("mu_upper", AxisArray([])),
        ("efficiency",AxisArray([])), ("e_set", AxisArray([]))]),
    transformers_t= Dict{String,AxisArray}( [("q1", AxisArray([])), ("q0", AxisArray([])), ("p0", AxisArray([])),
        ("p1", AxisArray([])), ("mu_upper", AxisArray([])), ("mu_lower", AxisArray([])),
        ("s_max_pu", AxisArray([]))]),
    snapshots=AxisArray([Any[]], [:t])

    )
    network_mutable(
        buses, generators, loads, lines, links, storage_units, stores, transformers, carriers,
        global_constraints,
        buses_t, generators_t, loads_t, lines_t, links_t, storage_units_t, stores_t, transformers_t,
        snapshots);
end
