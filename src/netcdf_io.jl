using DataFrames
using NetCDF

# This is not working yet. Rather experimental tries


function write_to_netcdffile(network, filename)
    tim=collect(1:size(network.snapshots)[1])
    components_t = [field for field=fieldnames(network) if String(field)[end-1:end]=="_t"]
    components_eq = [Symbol(String(field)[1:end-2]) for field=components_t]
    dyn_vs_stats = Dict(zip(components_t, components_eq))
    in(:loads_t,keys(dyn_vs_stats)) ? dyn_vs_stats[:loads_t]=:bus : nothing
    isfile("$filename") ? rm("$filename") : nothing
    for field_t in components_t
        field_eq = dyn_vs_stats[field_t]
        for df_name=keys(getfield(network,field_t))
            @show field_t, df_name
            if nrow(getfield(network,field_t)[df_name])>0
                df = getfield(network,field_t)[df_name][:,2:end]
                x_coords = collect(1:size(df)[1])
                nccreate("$filename","$field_t-$df_name", "time", tim, String(field_eq), x_coords)
                ncwrite(convert(Array{Float64}, df),"$filename","$field_t-$df_name")
            end
        end
    end
    components = [field for field=fieldnames(network) if String(field)[end-1:end]!="_t"]
    for field in components
        for attr in names(getfield(network, field))
            @show field, attr
            da = getfield(network,field)[attr]
            coords = collect(1:size(da)[1])
            nccreate("$filename","$field-$attr", String(field), coords)
            ncwrite(Array(da),"$filename","$field-$attr")
        end
    end
end
