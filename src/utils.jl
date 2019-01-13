using LightGraphs
using LightGraphs.SimpleGraphs
using Plots

include("compat.jl")

function time_dependent_components(network)
    fields = String.(fieldnames(network))
    components = []
    for field=fields
        if field[end-1:end] == "_t"
            push!(components, field)
        end
    end
    Symbol.(components)
end


function static_components(network)
    fields = String.(fieldnames(network))
    components = []
    for field=fields
        if field[end-1:end] != "_t"
            push!(components, field)
        end
    end
    Symbol.(components)
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


function align_component_order!(network)
    components_t = time_dependent_components(network)
    for comp=components_t
        order = Symbol.(getfield(network, Symbol(String(comp)[1:end-2]))[:name])
        for attr in keys(getfield(network, comp))
            if length(getfield(network,comp)[attr])==length(order)
                getfield(network,comp)[attr]= getfield(network,comp)[attr][:, order]
            end
        end
    end
end

idx(dataframe) = Dict(zip(dataframe[:name], Iterators.countfrom(1)))


rev_idx(dataframe) = Dict(zip(Iterators.countfrom(1), dataframe[:name]))


idx_by(dataframe, col, values) = select_by(dataframe, col, values)[:idx]


function select_by(dataframe, col, selector)
    if length(findin(dataframe[col], selector))==0
        return dataframe[repeat(Bool[false],outer=nrow(dataframe)) , :]
    else
        mdict = Dict(zip(dataframe[col], Iterators.countfrom(1)))
        ids = Array{Int,1}(0)
        for i in selector
            push!(ids, mdict[i])
        end
        dataframe[ids,:]
    end
end


select_names(a, b) = select_by(a, :name, b)


function append_idx_col!(dataframe)
    if typeof(dataframe)==Vector{DataFrames.DataFrame}
        for df in dataframe
            df[:idx] = collect(1:nrow(df))
        end
    else
        dataframe[:idx] = collect(1:nrow(dataframe))
    end
end


function get_switchable_as_dense(network, component, attribute, snapshots=0)
    snapshots==0 ? snapshots = network.snapshots : nothing
    T = nrow(snapshots)
    component_t = Symbol(component * "_t")
    component = Symbol(component)
    dense = DataFrame()
    if in(attribute, keys(getfield(network, component_t)))
        dense = getfield(network, component_t)[attribute]
    end
    cols = Symbol.(getfield(network, component)[:name])
    not_included = String.(setdiff(cols, names(dense)))
    if length(not_included)>0
        attribute = Symbol.(attribute)
        df = select_names(getfield(network, component), not_included)
        df = names!(DataFrame(repmat(transpose(Array(df[attribute])), T)),
                Symbol.(not_included))
        dense = [dense df]
    end
    return dense[cols]
end

"""Returns a function which selects time dependent variables either
from the series network.component_t or from static valuenetwork.component"""
function select_time_dep(network, component, attribute; components=0)

    component_t = Symbol(component * "_t")
    component = Symbol(component)
    attribute_s = Symbol(attribute)

    df = getfield(network, component_t)[attribute]

    if components==0
        c_names = Symbol.(getfield(network, component)[:,:name])
        s = getfield(network, component)[:,attribute_s]
    else
        c_names = Symbol.(getfield(network, component)[components,:name])
        s = getfield(network, component)[components,attribute_s]
    end

    df_names = names(df)

    return (t,i) -> in(c_names[i], df_names) ? df[t, c_names[i]] : s[i]
end


function calculate_dependent_values!(network)
    function set_default(dataframe, col, default)
        !in(col, names(dataframe)) ? dataframe[col] = default : nothing
    end

    #buses
    defaults = [(:v_nom, 1.)]
    for (col, default) in defaults
        set_default(network.buses, col, default)
    end

    # generators
    defaults = [(:p_nom_extendable, false), (:p_nom_max, Inf),(:commitable, false),
                (:p_min_pu, 0), (:p_max_pu, 1), (:p_nom_min, 0),(:capital_cost, 0),
                (:min_up_time, 0), (:min_down_time, 0), (:initial_status, true),
                (:p_nom, 0.),(:marginal_cost, 0),(:p_nom_opt, 0.)]
    for (col, default) in defaults
        set_default(network.generators, col, default)
    end

    # lines
    network.lines[:v_nom]=select_names(network.buses, network.lines[:bus0])[:v_nom]
    defaults = [(:s_nom_extendable, true), (:s_nom_min, 0),(:s_nom_max, Inf), (:s_nom, 0.),
                (:capital_cost, 0), (:g, 0), (:s_nom_ext_min, 0.1), (:s_max_pu, 1.0),] 
    for (col, default) in defaults
        set_default(network.lines, col, default)
    end
    
    network.lines[:x_pu] = network.lines[:x]./(network.lines[:v_nom].^2)
    network.lines[:r_pu] = network.lines[:r]./(network.lines[:v_nom].^2)
    #network.lines[:b_pu] = network.lines[:b].*network.lines[:v_nom].^2
    #network.lines[:g_pu] = network.lines[:g].*network.lines[:v_nom].^2

    # links
    defaults = [(:p_nom_extendable, false), (:p_nom_max, Inf), (:p_min_pu, 0),
                (:p_max_pu, 1),(:p_nom_min, 0), (:p_nom_max, Inf), (:capital_cost, 0),
                (:marginal_cost, 0), (:p_nom, 0.), (:efficiency, 1)]
    for (col, default) in defaults
        set_default(network.links, col, default)
    end

    # storage_units
    defaults = [(:p_nom_min, 0), (:p_nom_max, Inf), (:p_min_pu, -1),
                (:p_max_pu, 1), (:marginal_cost, 0), (:efficiency_store, 1),
                (:cyclic_state_of_charge, false),
                (:state_of_charge_initial, 0.), (:p_nom, 0.),
                (:efficiency_dispatch, 1), (:inflow, 0)]
    for (col, default) in defaults
        set_default(network.storage_units, col, default)
    end

    # stores
    defaults = [(:e_nom_min, 0), (:e_nom_max, Inf), (:e_min_pu, -1),
                    (:e_max_pu, 1), (:marginal_cost, 0), (:efficiency_store, 1),
                    (:efficiency_dispatch, 1),(:inflow, 0), (:e_nom, 0.)]
    for (col, default) in defaults
        set_default(network.stores, col, default)
    end

    # loads_t
    for df_name=keys(network.loads_t)
        if nrow(network.loads_t[df_name])>1
            for bus=[bus for bus in network.loads[:name] if
                !in(Symbol(bus), names(network.loads_t[df_name]))]
                set_default(network.loads_t[df_name], bus, 0)
            end
        end
    end

end

"""Converts the power system network into a Graph of LightGraphs.jl"""
function to_graph(network)
    busidx = idx(network.buses)
    g = DiGraph(length(busidx))
    for l in eachrow(network.lines)
        add_edge!(g, busidx[l[:bus0]], busidx[l[:bus1]] )
    end
    for l in eachrow(network.links)
        add_edge!(g, busidx[l[:bus0]], busidx[l[:bus1]] )
    end
    return g
end

"""Calculates the incidence matrix of a network"""
function incidence_matrix(network)
    busidx = idx(network.buses)
    lines = network.lines
    K = zeros(nrow(network.buses),nrow(lines))
    for l in 1:size(K)[2]
        K[busidx[lines[l,:bus0]],l] = 1
        K[busidx[lines[l,:bus1]],l] = -1
    end
    return K
end


"""Calculates the laplace matrix of a network"""
function laplace_matrix(network)
    K = incidence_matrix(network)
    return K*K'
end


"""Calculates the susceptance-weighted laplace matrix of a network"""
function weighted_laplace_matrix(network)
    K = incidence_matrix(network)
    B = susceptance_matrix(network)
    return K*B*K'
end


"""Returns the susceptance matrix of a network"""
function susceptance_matrix(network)
    return diagm(network.lines[:x].^(-1))
end


"""Calculates the ptdf matrix of a network"""
function ptdf_matrix(network)
    K = incidence_matrix(network)
    B = susceptance_matrix(network)
    L = weighted_laplace_matrix(network)
    H = B * K' * pinv(L)
    return H .- H[:,1]
end


"""Calculates the cycle basis of a network"""
function get_cycles(network)
    busidx = idx(network.buses)
    elist = [(busidx[l[:bus0]], busidx[l[:bus1]]) for l in eachrow(network.lines)]
    g = SimpleGraph(length(busidx))
    for e in elist add_edge!(g,e) end
    cycle_basis(g)
end


# TODO: not used
# following lift-and-project relaxation of Taylor2015a
function get_shortest_line_paths(network)
    busidx = idx(network.buses)

    # create elist only for lines with s_nom > 0
    elist = Tuple{Int64,Int64}[]
    for l in eachrow(network.lines)
        if l[:s_nom] > 0
            push!(elist, (busidx[l[:bus0]], busidx[l[:bus1]]))
        end
    end

    # create graph
    g = SimpleGraph(length(busidx))
    for e in elist add_edge!(g,e) end

    # set weights
    weights = zeros(nrow(network.lines), nrow(network.lines))
    for l in eachrow(network.lines)
        weights[busidx[l[:bus0]],busidx[l[:bus1]]] = l[:s_nom_step] * l[:x_step]
        weights[busidx[l[:bus1]],busidx[l[:bus0]]] = l[:s_nom_step] * l[:x_step]
    end

    line_paths = Array{Int64,1}[]
    for l in eachrow(network.lines)

        path_lg = a_star(g, busidx[l[:bus0]], busidx[l[:bus1]], weights)
        path = [(p.src,p.dst) for p in path_lg]
    
        lines_on_path = Int64[]
        for p in path
            for i=1:length(elist)
                if length(union(p, elist[i])) == length(p)
                    push!(lines_on_path, i)
                end
            end
        end

        push!(line_paths, lines_on_path)
    end

    return line_paths
end


function row_sum(df, row_id)
    if length(df[row_id,:]) == 0
        return 0.
    else
        return sum([df[row_id,i] for i in 1:length(df[row_id,:])])
    end
end


# TODO: generalise for multiple solvers; currently only works for Gurobi
"""Returns the irreducible infeasible set of constraints of an infeasible model"""
function get_iis(m::JuMP.Model)
    grb_model = MathProgBase.getrawsolver(internalmodel(m))
    num_constrs = Gurobi.num_constrs(grb_model)
    Gurobi.computeIIS(grb_model)
    iis_constrs = Gurobi.get_intattrarray(grb_model, "IISConstr",  1, num_constrs)
    m.linconstr[find(iis_constrs)]
end


# TODO: generalise for multiple solvers; currently only works for Gurobi
"""Returns the slack of each constraint of an optimised model"""
function get_slack(m::JuMP.Model)
    Gurobi.get_dblattrarray( m.internalModel.inner, "Slack", 1, Gurobi.num_constrs(m.internalModel.inner))
end


""""Returns the set of active constraints (zero slack) of an optimised model"""
function get_active_constraints(m::JuMP.Model)
    
    slack = get_slack(m)
    
    function within_tolerance(x)
        (x <= 1e-9 && x >= -1e-9)
    end
    
    sort(collect(Set(find(within_tolerance,slack))))
end


"""Prints the set of active constraints with slack (should be zero) to the console"""
function print_active_constraints!(m::JuMP.Model)

    slack = get_slack(m)
    active_constraints = get_active_constraints(m)

    for ac in active_constraints
        println(round(slack[ac],5),"\t\t",m.linconstr[ac])
    end

end


""""Returns the set of inactive constraints (non-zero slack) of an optimised model"""
function get_inactive_constraints(m::JuMP.Model; slack_filter=1e-9)
    
    slack = get_slack(m)
    
    function within_filter(x)
        (x >= slack_filter || x <= -slack_filter)
    end
    
    sort(collect(Set(find(within_filter,slack))))
end


"""Prints the set of inactive constraints with slack (should be non-zero) to the console"""
function print_inactive_constraints!(m::JuMP.Model; slack_filter=1e-9)

    slack = get_slack(m)
    inactive_constraints = get_inactive_constraints(m, slack_filter=slack_filter)

    for iac in inactive_constraints
        println(round(slack[iac],5),"\t\t",m.linconstr[iac])
    end

end


# simple approach; preferred approach is preprocessing the snapshots e.g. with the tsam package
"""Samples snapshots at specified rate and adapts snapshot weights uniformly"""
function set_snapshots_sampling!(network, sampling_rate)
    T = nrow(network.snapshots)
    rows_to_delete = setdiff(collect(1:T),collect(1:sampling_rate:T))
    network.snapshots = deleterows!(network.snapshots,rows_to_delete)
    network.generators_t["p_max_pu"] = deleterows!(network.generators_t["p_max_pu"],rows_to_delete)
    network.loads_t["p_set"] = deleterows!(network.loads_t["p_set"],rows_to_delete)
    network.snapshots[:weightings] *= sampling_rate
end


# TODO: support for multiple line types; currently frequently used line type is hard coded
"""Sets the maximum number of extendable circuits for each line"""
function set_maximum_extendable_circuits!(network; additional_num_parallel=0, extension_factor=1)
    for i=1:nrow(network.lines)
        # allow some buffer to avoid rounding issues
        network.lines[:s_nom_max][i] =  extension_factor*network.lines[:s_nom][i] + additional_num_parallel*1698.11
    end
end


"""Returns the set of candidate lines (number of additional circuits) for each line"""
function line_extensions_candidates(network)
    # candidates set always includes 0
    candidates = Array{Int64,1}[]
    lines = network.lines
    N_ext_LN = sum(.!(.!lines[:s_nom_extendable]))
    for l=1:N_ext_LN
        if lines[:s_nom_max][l] != Inf
            max_extension_float = (lines[:s_nom_max][l] / lines[:s_nom][l] - 1) * lines[:num_parallel][l]
            max_extension = floor(max_extension_float) 
            push!(candidates,[i for i=0:max_extension]) # starts from 0 as no extension is also an
        else
            # fallback option
            push!(candidates,[i for i=0:2]) # starts from 0 as no extension is also an
        end
    end
    return candidates
end


"""Returns the constraint matrix of a model.
If nz=true: entries have value 0 or 1
If nz=false: entries have value of corresponding parameter
"""
function constraintmatrix(model::JuMP.Model; nz::Bool=false)
    constraintmatrix = MathProgBase.getconstrmatrix(internalmodel(model))
    cm = deepcopy(constraintmatrix)
    if nz
        nzcm = findnz(cm)
        for (i, j) in zip(nzcm[1],nzcm[2])
            cm[i,j] = 1.0
        end
    end
    return cm
end

"""Plots the constraint matrix"""
function plot_cm_nonzeroentries(model::JuMP.Model; cm=nothing)
    if cm == nothing
        cm = constraintmatrix(model, nz=true)
        colSwitch = colswitch_cm(model)
        println(colSwitch)
        cm = full(cm)[:,colSwitch]
    end
    xs = [string("x",i) for i = 1:size(cm)[1]]
    ys = [string("y",i) for i = 1:size(cm)[2]]
    heatmap(xs,ys,cm',aspect_ratio=1, color=:Greys, size=(1300,800), title="Non-zero entries of constraint matrix")
end


"""Plots the value distribution of the constraint matrix"""
function plot_cm_valuedistribution(model::JuMP.Model; cutoff=1e6, cm=nothing)
    if cm == nothing
        cm = constraintmatrix(model, nz=false)
    end
    elements = findnz(cm)[3]
    min=minimum(abs.(elements))
    max=maximum(abs.(elements))
    filter!(x -> x <= cutoff, elements)
    filter!(x -> x >= -cutoff, elements)
    histogram(elements,nbins=100, title="Distribution of values in constraint matrix -- absolute value range: [$min, $max]",
        xlabel="value / coefficient", ylabel="frequency", size=(1400,800), legend=false, color=:grays)
end


"""Switches columns of constraint matrix for organised plotting"""
function colswitch_cm(model::JuMP.Model)
    model.objDict # produces model.colNames
    variables = model.colNames
    println(model.colNames)
    switches = []
    T = size(model[:LN_ext])[2]
    push!(switches,find(x->x==true, .![contains(i, ",") for i in variables]))
    for t=1:T
        push!(switches,find(x->x==true, [contains(i, ",$t]") for i in variables]))
    end
    return reverse(vcat(switches...))
end


"""Calculates the Big-M parameters for each line using the assumed maximum angle difference."""
function bigm(cnstr::Symbol, network; max_angle_diff::Float64=pi/6)
                
    lines = network.lines
    candidates = line_extensions_candidates(network)
    init_c = lines[:num_parallel]
    init_x = lines[:x_pu]
    init_s_nom = lines[:s_nom]
    max_c = maximum.(candidates)

    if cnstr == :flows_upper
        extreme_c = minimum.(candidates)
    elseif cnstr == :flows_lower
        extreme_c = maximum.(candidates)
    else
        error("Function bigm(.) not defined for constraint $constr.")
    end
    
    bigm = ( extreme_c ./ init_c .+ 1) .* ( max_angle_diff ./ init_x ) 
        .+ ( max_c ./ init_c .+ 1 ) .* init_s_nom 
    
    length(bigm) != nrow(lines) ? error("Sizes of Big-M parameter and number of lines do not match!") : nothing
    
    return 10*bigm

end


"""Retrieves all variable symbols (names) of a model"""
function getvariables(m::JuMP.Model)
    return [k for (k,v) in m.objDict 
             if issubtype(eltype(v), JuMP.Variable) &&
             k!= :ALPHA &&
             size(m[k])[1]!=0
           ]
end


"""Retrieves all constraint symbols (names) of a model"""
function getconstraints(m::JuMP.Model)
    objs = [k for (k,v) in m.objDict if issubtype(eltype(v), JuMP.ConstraintRef)]
    constraints = []
    for obj in objs
        empty = false
        try
            m[obj][1,1]
        catch
            empty = true
        end
        if !empty
            push!(constraints, obj)
        end
    end
    return constraints
end


JuMP.rhs(constraint::JuMP.ConstraintRef) = JuMP.rhs(LinearConstraint(constraint))


JuMP.rhs(constraint::Array) = reshape([JuMP.rhs(constraint[1])],1,1) # special case for global constraints


JuMP.rhs(constraints::JuMP.JuMPArray{JuMP.ConstraintRef}) = JuMP.JuMPArray(JuMP.rhs.(constraints.innerArray), constraints.indexsets)


"""Computes the constant term of a Benders cut using the uncoupled constraints of the subproblems."""
function get_benderscut_constant(m::JuMP.Model, uncoupled_constraints::Array{Any,1})
    constant = 0
    for constr in uncoupled_constraints
        try
            constant += dot(getdual(m[constr])[:,:],JuMP.rhs(m[constr])[:,:])
        catch
            constant += dot(getdual(m[constr])[:],JuMP.rhs(m[constr])[:])
        end
    end
    return constant
end

"""Computes the constant term of a Benders cut using the uncoupled constraints of the subproblems."""
function get_benderscut_constant(m::Array{JuMP.Model,1}, uncoupled_constraints::Array{Any,1})
    sum = 0
    for i=1:length(m)
        sum += get_benderscut_constant(m[i], uncoupled_constraints)
    end
    return sum
end


function getduals(m::Array{JuMP.Model,1}, cnstr::Symbol; filter_b::Bool=false)
    z = vcat(getfield.(getdual.(getindex.(m,cnstr)), :innerArray)...)
    filter_b ? z[(z.<1e-3).&(z.>-1e-3)] = 0.0 : nothing
    return z
end


function getduals_flows(m::Array{JuMP.Model,1}, cnstr::Symbol; filter_b::Bool=false)
    x=getdual(getindex.(m, cnstr))
    if filter_b
        for t=1:length(x)
             for (l,c,ts) in keys(x[t])
                z = x[t][l,c,ts]
                if z>-1e-4 && z<1e-4
                    x[t][l,c,ts] = 0.0
                end
             end
        end
    end
    return x
end


"""Specifies the rescaling factors applied to different equations of the problem"""
function rescaling_factors(rescaling::Bool)
    # TODO: adapt rescaling factors
    dict = Dict(
        :approx_restarget => 1e-3,
        :bounds_G => 1e4,
        :bounds_LN => 1e3,
        :bounds_LK => 1e2,
        :flows => 1,
        :benderscut => 1e-3,
        :objective => 1
    )
    !rescaling ? for k in keys(dict) dict[k] = 1 end : nothing
    return dict
end


"""Filters extremely small time series values (e.g. p_max_pu) below a threshold"""
function filter_timedependent_extremes!(z, threshold::Float64)
    for c in z.df.columns[2:end] 
        c[(c.<threshold).&(c.!=0)] = threshold 
    end
    return z.df
end


"""Organises optimisation output"""
function write_optimalsolution(network, m::JuMP.Model; sm=nothing, joint::Bool=true)

    joint ? sm = m : nothing

    if (typeof(sm) != Array{JuMP.Model,1})
        sm = [sm]
    end

    # shortcuts
    buses = network.buses
    lines = network.lines
    generators = network.generators
    links = network.links
    stores = network.stores
    storage_units = network.storage_units

    N = nrow(network.buses)
    L = nrow(network.lines)
    T = nrow(network.snapshots)

    fix_gens_b = (.!generators[:p_nom_extendable])
    ext_gens_b = .!fix_gens_b
    fix_lines_b = (.!lines[:s_nom_extendable])
    ext_lines_b = .!fix_lines_b
    fix_links_b = .!links[:p_nom_extendable]
    ext_links_b = .!fix_links_b
    fix_sus_b = .!storage_units[:p_nom_extendable]
    ext_sus_b = .!fix_sus_b
    fix_stores_b = .!stores[:e_nom_extendable]
    ext_stores_b = .!fix_stores_b

    G = [hcat(getindex.(sm, :G_fix)...); hcat(getindex.(sm, :G_ext)...)]
    LN = [hcat(getindex.(sm, :LN_fix)...); hcat(getindex.(sm, :LN_ext)...)]
    LK = [hcat(getindex.(sm, :LK_fix)...); hcat(getindex.(sm, :LK_ext)...)]
    SU_dispatch = [hcat(getindex.(sm, :SU_dispatch_fix)...); hcat(getindex.(sm, :SU_dispatch_ext)...)]
    SU_store = [hcat(getindex.(sm, :SU_store_fix)...); hcat(getindex.(sm, :SU_store_ext)...)]
    SU_soc = [hcat(getindex.(sm, :SU_soc_fix)...); hcat(getindex.(sm, :SU_soc_ext)...)]
    SU_spill = [hcat(getindex.(sm, :SU_spill_fix)...); hcat(getindex.(sm, :SU_spill_ext)...)]
    ST_dispatch = [hcat(getindex.(sm, :ST_dispatch_fix)...); hcat(getindex.(sm, :ST_dispatch_ext)...)]
    ST_store = [hcat(getindex.(sm, :ST_store_fix)...); hcat(getindex.(sm, :ST_store_ext)...)]
    ST_soc = [hcat(getindex.(sm, :ST_soc_fix)...); hcat(getindex.(sm, :ST_soc_ext)...)]
    ST_spill = [hcat(getindex.(sm, :ST_spill_fix)...); hcat(getindex.(sm, :ST_spill_ext)...)]

    capex_gep_cost = dot(generators[ext_gens_b,:capital_cost], getvalue(m[:G_p_nom]))
    + dot(generators[fix_gens_b,:capital_cost], generators[fix_gens_b,:p_nom])

    capex_tep_cost = dot(lines[ext_lines_b,:capital_cost], getvalue(m[:LN_s_nom])) 
    + dot(lines[fix_lines_b,:capital_cost], lines[fix_lines_b,:s_nom]) 

    lines = [lines[fix_lines_b,:]; lines[ext_lines_b,:]]
    orig_gen_order = network.generators[:name]
    generators = [generators[fix_gens_b,:]; generators[ext_gens_b,:] ]
    links = [links[fix_links_b,:]; links[ext_links_b,:]]
    storage_units = [storage_units[fix_sus_b,:]; storage_units[ext_sus_b,:]]
    stores = [stores[fix_stores_b,:]; stores[ext_stores_b,:]]

    opex_cost = sum(network.snapshots[:weightings][t]*dot(generators[:marginal_cost], getvalue(G[:,t])) for t=1:T)

    # tep costs
    tep_cost = dot(lines[ext_lines_b,:capital_cost], getvalue(m[:LN_s_nom]))

    # write results
    generators[:p_nom_opt] = deepcopy(generators[:p_nom])
    fix_gens_b_reordered = (.!generators[:p_nom_extendable])
    ext_gens_b_reordered = .!fix_gens_b_reordered
    generators[ext_gens_b_reordered,:p_nom_opt] = getvalue(m[:G_p_nom])
    network.generators = generators
    network.generators_t["p"] = names!(DataFrame(transpose(getvalue(G))), Symbol.(generators[:name]))
    network.generators = select_names(network.generators, orig_gen_order)

    orig_line_order = network.lines[:name]
    network.lines = lines
    lines[:s_nom_opt] = deepcopy(lines[:s_nom])
    network.lines[ext_lines_b,:s_nom_opt] = getvalue(m[:LN_s_nom])
    network.lines_t["p0"] = names!(DataFrame(transpose(getvalue(LN))), Symbol.(lines[:name]))
    network.lines = select_names(network.lines, orig_line_order)

    if nrow(links)>0
        orig_link_order = network.links[:name]
        network.links = links
        links[:p_nom_opt] = deepcopy(links[:p_nom])
        network.links[ext_links_b,:p_nom_opt] = getvalue(m[:LK_p_nom])
        network.links_t["p0"] = names!(DataFrame(transpose(getvalue(LK))), Symbol.(links[:name]))
        network.links = select_names(network.links, orig_link_order)

    end
    if nrow(storage_units)>0
        orig_sus_order = network.storage_units[:name]
        network.storage_units = storage_units

        storage_units[:p_nom_opt] = deepcopy(storage_units[:p_nom])
        network.storage_units[ext_sus_b,:p_nom_opt] = getvalue(m[:SU_p_nom])
        network.storage_units_t["spill"] = names!(DataFrame(transpose(getvalue(SU_spill))),
                            Symbol.(storage_units[:name]))
        network.storage_units_t["p"] = names!(DataFrame(transpose(getvalue(SU_dispatch .- SU_store))),
                            Symbol.(storage_units[:name]))
        network.storage_units_t["state_of_charge"] = names!(DataFrame(transpose(getvalue(SU_soc))),
                            Symbol.(storage_units[:name]))
        network.storage_units = select_names(network.storage_units, orig_sus_order)
    end

    align_component_order!(network)

    # print total system code
    println("\nObjective:\t$(m.objVal)")
    println("TEP share:\t$(tep_cost/m.objVal*100) %")

    generators = [generators[fix_gens_b_reordered,:]; generators[ext_gens_b_reordered,:] ]
    nonnull_carriers = network.carriers[network.carriers[:co2_emissions].!=0, :][:name]
    carrier_index(carrier) = findin(generators[:carrier], carrier)

    # renewable energy target
    null_carriers = network.carriers[network.carriers[:co2_emissions].==0,:][:name]
    res = sum(sum(network.snapshots[:weightings][t]*getvalue(G[carrier_index(null_carriers),t]) for t=1:T)) / 
    sum(sum(network.snapshots[:weightings][t]*getvalue(G[:,t]) for t=1:T))

    println("RES share:\t$(res*100) %")

    # curtailment
    sum_of_dispatch = sum(network.snapshots[:weightings][t]*getvalue(G[findin(generators[:name],string.(network.generators_t["p_max_pu"].colindex.names)),t]) for t=1:T)
    p_nom_opt = generators[:p_nom_opt][findin(generators[:name], string.(network.generators_t["p_max_pu"].colindex.names))]
    ren_gens_i = findin(network.generators[:name], string.(network.generators_t["p_max_pu"].colindex.names))
    ren_gens_b = [in(i,ren_gens_i) ? true : false for i=1:size(network.generators)[1]]
    fix_ren_gens_b = .!network.generators[:p_nom_extendable][ren_gens_b]
    p_max_pu = network.generators_t["p_max_pu"][:,2:end]
    p_max_pu = [p_max_pu[:,fix_ren_gens_b] p_max_pu[:,.!fix_ren_gens_b]]
    p_max_pu = convert(Array,p_max_pu)
    sum_of_p_max_pu = sum(network.snapshots[:weightings][t]*p_max_pu[t,:] for t=1:T)
    curtailment = 1 - ( sum(sum_of_dispatch) / dot(p_nom_opt,sum_of_p_max_pu) )

    println("Curtailment:\t$(curtailment*100) %")

    println("")

end