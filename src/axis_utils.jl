
idx(dataframe) = Dict(zip(dataframe.axes[1].val, Iterators.countfrom(1)))
rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe.axes[1].val))

zsum(array) = length(array)>0 ? sum(array) : 0.
zdot(v1,v2) = length(v1)>0 ? dot(v1,v2) : 0.


"""
Use this funtion to assign a new column or row with given values.
Add a new row with keywordargument axis=1, for column set axis=2.
"""
function assign(df::AxisArray, values, index; axis=2)
    # This could also be done with AxisArrays.merge which breaks however with type 
    # Any. 
    # check if index is already in manipulated axis of df
    in(typeof(index),[String, Symbol, DateTime]) ? index = [index] : nothing
    if index[1] ∈ df.axes[axis].val
        df = deepcopy(df)
        axis == 1 ? df[index ,:] .= values' : nothing
        axis == 2 ? df[:, index] .= values : nothing
        df
    else
        if axis==1
            values = isa(values, Array) ? values : fill(values, (size(df)[2])) 
            AxisArray([df.data; values'], 
                      Axis{axisnames(df)[1]}(append!(copy(axes(df)[1].val), index)),  
                      axes(df)[2])
        elseif axis == 2
            values = isa(values, Array) ? values : fill(values, (size(df)[1])) 
            AxisArray([df.data values], 
                        axes(df)[1], 
                        Axis{axisnames(df)[2]}(append!(copy(axes(df)[2].val), index))) 
        end
    end
end

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

