

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

