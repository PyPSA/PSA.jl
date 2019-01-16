module PSA

using DataFrames, CSV, AxisArrays, NCDatasets, NamedTuples, Dates

export Network, import_nc, export_nc

include("lopf_pathway.jl")
# include("lopf.jl")


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
    buses_t::NamedTuple
    generators_t::NamedTuple
    loads_t::NamedTuple
    lines_t::NamedTuple
    links_t::NamedTuple
    storage_units_t::NamedTuple
    stores_t::NamedTuple
    transformers_t::NamedTuple
    snapshots::Array
    snapshot_weightings::AxisArray
    name::String
end

function Network(
# static
    buses = AxisArray(repeat(Any[], inner=(1,16)), 
            Axis{:row}(String[]),
            string.([:name, :v_nom, :type, :x, :y, :carrier, :v_mag_pu_set, :v_mag_pu_min,
            :v_mag_pu_max, :control, :sub_network, :p, :q,
            :v_mag_pu, :v_ang, :marginal_price])),

    carriers = AxisArray(repeat(Any[], inner=(1,2)), 
                Axis{:row}(String[]),
                string.([:name, :co2_emissions])),

    generators = AxisArray(repeat(Any[], inner=(1,31)), 
            Axis{:row}(String[]),
            string.([:name, :bus, :control, :type, :p_nom, :p_nom_extendable, :p_nom_opt,
            :p_nom_min, :p_nom_max, :p_min_pu, :p_max_pu, :p_set, :q_set,
            :sign, :carrier, :marginal_cost, :capital_cost, :efficiency,
            :committable, :start_up_cost, :shut_down_cost, :min_up_time,
            :min_down_time, :initial_status, :ramp_limit_up, :ramp_limit_down,
            :ramp_limit_start_up, :ramp_limit_shut_down, :p, :q, :status])),

    global_constraints = AxisArray(repeat(Any[], inner=(1,6)), 
            Axis{:row}(String[]),
            string.([:name, :type, :carrier_attribute, :sense, :constant, :mu])),

    line_types = AxisArray(repeat(Any[], inner=(1,9)), 
            Axis{:row}(String[]),
            string.([:name, :f_nom, :r_per_length, :x_per_length, :c_per_length,
            :i_nom, :mounting, :cross_section, :references])),

    lines = AxisArray(repeat(Any[], inner=(1,33)), 
            Axis{:row}(String[]),
            string.([:name, :bus0, :bus1, :type, :x, :r, :g, :b,
            :s_nom, :s_nom_extendable, :s_nom_opt, :s_nom_min, :s_nom_max,
            :s_max_pu, :capital_cost, :length, :terrain_factor, :num_parallel,
            :v_ang_min, :v_ang_max, :sub_network, :p0, :q0, :p1, :q1,
            :x_pu, :r_pu, :g_pu, :b_pu, :x_pu_eff, :r_pu_eff,
            :mu_lower, :mu_upper])),

    links = AxisArray(repeat(Any[], inner=(1,21)), 
            Axis{:row}(String[]),
            string.([:name, :bus0, :bus1, :type, :efficiency, :p_nom, :p_nom_extendable,
            :p_nom_min, :p_nom_max, :p_set, :p_min_pu, :p_max_pu, :capital_cost,
            :marginal_cost, :length, :terrain_factor, :p0, :p1,
            :p_nom_opt, :mu_lower, :mu_upper])),

    loads = AxisArray(repeat(Any[], inner=(1,8)), 
            Axis{:row}(String[]),
            string.([:name, :bus, :type, :p_set, :q_set, :sign, :p, :q])),

    shunt_impendances = AxisArray(repeat(Any[], inner=(1,9)), 
            Axis{:row}(String[]),
            string.([:name, :bus, :g, :b, :sign, :p, :q, :g_pu, :b_pu])),

    storage_units = AxisArray(repeat(Any[], inner=(1,29)), 
            Axis{:row}(String[]),
            string.([:name, :bus, :control, :type, :p_nom, :p_nom_extendable, :p_nom_min,
            :p_nom_max, :p_min_pu, :p_max_pu, :p_set, :q_set, :sign, :carrier,
            :marginal_cost, :capital_cost, :state_of_charge_initial, :state_of_charge_set,
            :cyclic_state_of_charge, :max_hours, :efficiency_store, :efficiency_dispatch,
            :standing_loss, :inflow, :p, :q, :state_of_charge, :spill, :p_nom_opt])),

    stores = AxisArray(repeat(Any[], inner=(1,22)), 
            Axis{:row}(String[]),
            string.([:name, :bus, :type, :e_nom, :e_nom_extendable, :e_nom_min,
            :e_nom_max, :e_min_pu, :e_max_pu, :e_initial, :e_cyclic, :p_set,
            :cyclic_state_of_charge, :q_set, :sign, :marginal_cost,
            :capital_cost, :standing_loss, :p, :q, :e, :e_nom_opt])),

    transformer_types = AxisArray(repeat(Any[], inner=(1,16)), 
            Axis{:row}(String[]),
            string.([:name, :f_nom, :s_nom, :v_nom_0, :v_nom_1, :vsc, :vscr, :pfe,
            :i0, :phase_shift, :tap_side, :tap_neutral, :tap_min, :tap_max, :tap_step, 
            :references])),

    transformers = AxisArray(repeat(Any[], inner=(1,36)), 
            Axis{:row}(String[]),
            string.([:name, :bus0, :bus1, :type, :model, :x, :r, :g, :b, :s_nom,
            :s_nom_extendable, :s_nom_min, :s_nom_max, :s_max_pu, :capital_cost,
            :num_parallel, :tap_ratio, :tap_side, :tap_position, :phase_shift,
            :v_ang_min, :v_ang_max, :sub_network, :p0, :q0, :p1, :q1, :x_pu,
            :r_pu, :g_pu, :b_pu, :x_pu_eff, :r_pu_eff, :s_nom_opt, :mu_lower, :mu_upper])),

# time_dependent
    buses_t = (marginal_price = AxisArray([]), v_ang =  AxisArray([]),
            v_mag_pu_set = AxisArray([]), q =  AxisArray([]),
            v_mag_pu =  AxisArray([]), p =  AxisArray([]), 
            p_max_pu =  AxisArray([])),
    generators_t = (
            marginal_price = AxisArray([]), v_ang =  AxisArray([]),
            v_mag_pu_set = AxisArray([]), q =  AxisArray([]),
            v_mag_pu =  AxisArray([]), p =  AxisArray([]),
            p_min_pu =  AxisArray([]), p_max_pu =  AxisArray([]),
            p_nom_opt =  AxisArray([])),
    loads_t = (q_set = AxisArray([]), p_set =  AxisArray([]),
            q =  AxisArray([]), p =  AxisArray([])),
    lines_t = (q0 = AxisArray([]), q1 =  AxisArray([]),
            p0 = AxisArray([]), p1 =  AxisArray([]),
            mu_lower =  AxisArray([]), mu_upper =  AxisArray([]),
            s_nom_opt =  AxisArray([])),
    links_t = (p_min_pu = AxisArray([]), p_max_pu =  AxisArray([]),
            p0 = AxisArray([]), p1 =  AxisArray([]),
            mu_lower =  AxisArray([]), mu_upper =  AxisArray([]),
            efficiency = AxisArray([]), p_set =  AxisArray([]),
            p_nom_opt =  AxisArray([])),
    storage_units_t = (p_min_pu = AxisArray([]), p_max_pu =  AxisArray([]),
            inflow = AxisArray([]), mu_lower =  AxisArray([]), mu_upper =  AxisArray([]),
            efficiency = AxisArray([]), p_set =  AxisArray([])),
    stores_t=(e_min_pu = AxisArray([]), e_max_pu =  AxisArray([]),
            inflow = AxisArray([]), mu_lower =  AxisArray([]), mu_upper =  AxisArray([]),
            efficiency = AxisArray([]), e_set =  AxisArray([])),
    transformers_t= (q1 =  AxisArray([]), q0 =  AxisArray([]), p0 =  AxisArray([]),
            p1 =  AxisArray([]), mu_upper =  AxisArray([]), mu_lower =  AxisArray([]),
            s_max_pu =  AxisArray([])),
    snapshots=Array([]), 
    snapshot_weightings = AxisArray(Any[], Axis{:row}(String[])),
    name = "PSA network"
    )
    network_mutable(
        buses, generators, loads, lines, links, storage_units, stores, transformers, carriers,
        global_constraints,
        buses_t, generators_t, loads_t, lines_t, links_t, storage_units_t, stores_t, transformers_t,
        snapshots, snapshot_weightings, name);
end



# auxiliary for import_nc
function reformat(data::Array)
    if typeof(data) == Array{Union{Missing, Int8},1}
        Array{Bool}(data)
    elseif typeof(data) == Array{Union{Missing, Int64},1}
        Array{Float64}(data)
    else
        data    
    end
end

function import_nc(path)
    n = Network()
    ds = Dataset(path, )
    ds_keys = keys(ds) 
    stats = static_components(n) 
    push!(stats, "snapshots_weightings")
    found = ""
    for comp = string.(stats)
        if any(contains.(ds_keys, comp))
            found = found * "$comp, "
            if comp == "snapshots"
                data = ds["snapshots"][:]
                setfield!(n, Symbol(comp), Array(data))
            elseif comp == "snapshots_weightings"
                data = ds["snapshots_weightings"][:]
                setfield!(n, Symbol("snapshot_weightings"),
                AxisArray(data, Axis{:time}(ds["snapshots"][:])))
            else
                index = reformat(ds[comp * "_i"][:])
                props = ds_keys[findall(contains.(ds_keys,
                            Regex("$(comp)_(?!(i\$|t_))"))) ]
                cols = props
                setfield!(n, Symbol(comp),
                        AxisArray(  hcat([reformat(ds[key][:]) for key = props]...) ,
                        Axis{:row}(index), Axis{:col}(replace.(props, comp * "_", ""))) )
            end
        end
    end
    @info("The following static components were imported: \n $found")
    dyns = dynamic_components(n)
    found = ""
    for comp=dyns
        comp_stat = Symbol(String(comp)[1:end-2])
        for attr in keys(getfield(n, comp))
            if in("$(comp)_$attr", ds_keys)
                found = found * "$(comp)_$attr, "
                val = AxisArray( ds["$(comp)_$attr"][:]',
                        Axis{:snapshots}(n.snapshots), 
                        Axis{comp_stat}(ds["$(comp)_$(attr)_i"][:]))
                replace_attribute!(n, comp, attr, val)
            end
        end
    end
    n.name = ds.attrib["network_name"]
    close(ds)
    @info("The following dynamic components were imported: \n $(found)")
    n;
end

function catch_type(array)
    assert(ndims(array)==1)
    i = 0
    t = Missing
    while t == Missing
        i += 1
        t = typeof(array[i])
    end
    t==Bool ? Int8 : t
end


function export_nc(n, path)
    ispath(path) ? rm(path) : nothing
    ds = Dataset(path,"c")
    found = ""
    for comp = static_components(n)
        data = getfield(n, comp)
        comp = String(comp)
        if length(data)  > 0 
            found = found * "$comp, "
            if comp == "snapshots"
                # Dimension
                defDim(ds, comp, length(data))
                # Variable
                sn = defVar(ds, comp, Int64, (comp,))
                sn.attrib["units"] = "hours since $(data[1])"
                sn.attrib["calendar"] = "gregorian"
                sn[:] = 1:length(data)
            elseif comp == "snapshots_weightings"
                continue
            else
                # Dimension
                rows = data.axes[1].val
                string_length = maximum(length.(rows))
                ds.dim["$(comp)_i"] = size(rows)[1]
                # ds.dim["string$string_length"] = string_length

                # Variable
                index = defVar(ds, "$(comp)_i", String, ("$(comp)_i", ))
                index[:] = rows
                for col in data.axes[2].val
                    vartype = catch_type(data[:,col])
                    var = defVar(ds, "$(comp)_$col", vartype, ("$(comp)_i",))
                    vartype == Int8 ? var.attrib["dtype"] = "bool" : nothing 
                    vartype != String ? var[:] = collect(Missings.replace(data[:, col].data, NaN )) : nothing
                    vartype == String ? var[:] = data[:, col].data : nothing
                    var.attrib["FillValue"] = NaN
                end
            end
        end
    end
    @info("The following static components were exported: \n $found")
    found = ""
    for comp=dynamic_components(n)
        for attr in keys(getfield(n, comp))
            data = getfield(n, comp)[attr]
            # comp = string(comp)
            if length(data)  > 0 
                found = found * "$(string(comp)[1:end-2]), "
                # Dimension
                cols = data.axes[2].val
                string_length = maximum(length.(cols))
                ds.dim["$(comp)_$(attr)_i"] = length(cols)
                # ds.dim["string$string_length"] = string_length


                index = defVar(ds, "$(comp)_$(attr)_i", String,  
                        ("$(comp)_$(attr)_i",) )
                index[:] = cols
                var = defVar(ds, "$(comp)_$attr", Float64, 
                        ("$(comp)_$(attr)_i","snapshots",  ))
                var[:,:] = (data.data)'
            end
        end
    end
    #ds.attrib["_NCProperties"] = "..."
    ds.attrib["network_srid"] = ""
    ds.attrib["network_pypsa_version"] = "0.13.1"
    ds.attrib["network_name"] = n.name

    @info("The following dynamic components were exported: \n $(found)")
    close(ds)
end

# ------------------------------------------------------------------------------------------------

# auxiliary for import_csv
function df2axarray(df; indexname=:row, colname=:col, index=nothing, dtype=Any)
    # take first column as index 
    index == nothing ? index = string.(df[1]) : nothing
    cols = string.(names(df)[2:end])
    data = Array{dtype}(df[:,2:end])
    AxisArray(data, Axis{indexname}(index), Axis{colname}(cols))
end

# auxiliary for import_csv
to_datetime(stringarray) = DateTime.(stringarray, "y-m-d H:M:S")

function import_csv(folder)
    # reeeally slow right now because of the slow conversion from DataFrames.DataFrame to Array.
    # Especially for the dynamic data.
    # this would be better if one could set the import type in readdlm to only floats
    # or skip the first column which is the snapshots-string-column
    n = Network()
    !ispath("$folder") ? error("Path not existent") : nothing
    stats = static_components(n)
    for comp=stats
        if ispath("$folder/$comp.csv")
            if comp == :snapshots 
                sns = CSV.read("$folder/$comp.csv")
                sns = to_datetime(sns[:name])
                setfield!(n, comp, sns)
                continue
            end 
            df = CSV.read("$folder/$comp.csv"; truestring="True", falsestring="False")
            setfield!(n, Symbol(comp), df2axarray(df))
        end
    end
    dyns = dynamic_components(n)
    for comp_t=dyns
        for attr in keys(getfield(n, comp_t))
            comp = Symbol(String(comp_t)[1:end-2])
            if ispath("$folder/$comp-$attr.csv")
                getfield(n,comp_t)[attr] = (
                    df2axarray(CSV.read("$folder/$comp-$attr.csv"; 
                                truestring="True", falsestring="False"), 
                            index = n.snapshots, dtype=Float64) )
            end
        end
    end
    initializer = Network()
    for field=setdiff(fieldnames(n), fieldnames(initializer))
        setfield!(n, field, getfield(initializer, field))
    end
    n
end


function export_csv(n, folder)
    !ispath(folder) ? mkdir(folder) : nothing
    dyns = [field for field=fieldnames(n) if String(field)[end-1:end]=="_t"]
    for field_t in dyns
        for df_name=keys(getfield(n,field_t))
            if nrow(getfield(n,field_t)[df_name])>0
                field = Symbol(String(field_t)[1:end-2])
                writetable("$folder/$field-$df_name.csv", getfield(n,field_t)[df_name])
            end
        end
    end
    stats = [field for field=fieldnames(n) if String(field)[end-1:end]!="_t"]
    for field in stats
        writetable("$folder/$field.csv", getfield(n, field))
    end
end





end
