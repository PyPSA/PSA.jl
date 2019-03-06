using JuMP
using MathProgBase
using PowerModels
pm = PowerModels

include("utils.jl")

"""Build a linear optimal power flow model"""
function build_lopf(network, solver; rescaling::Bool=false,formulation="angles_linear",#::Union{String, Type{GenericPowerModel{F}}}="angles_linear",
                    investment_type::String="continuous",
                    blockmodel::Bool=false, benders::String="", snapshot_number=0, N_cuts=1, blockstructure::Bool=false) #where F <: pm.AbstractPowerFormulation

    # ---------------- #               
    # Function Outline #
    # ---------------- #

    # 1.        Initialize model and auxiliaries
    # 2.        Get base model from PowerModels.jl
    # 3.        Add all variables
    # 4. - 10.  Add generators,  lines, links, storage_units, shunts, stores to the model:
    #               .1 separate different types from each other
    #               .2 define number of different types
    #               .3 add variables to the model
    #               .4 set contraints for extendables
    #               .5 set charging constraints (storage_units and stores)
    # 11.       give power flow formulation
    # 12.       set global constraints
    # 13.       give objective function
    
    # ------------ #
    # Conventions: #
    # ------------ #
    
    #   - all variable names start with capital letters
    

# --------------------------------------------------------------------------------------------------------
# 1. Initialize model and auxiliaries
# --------------------------------------------------------------------------------------------------------

    # sanity checks
    snapshot_number>0 && benders!="slave" ? error("Can only specify one single snapshot for slave-subproblem!") : nothing
    blockmodel && benders!="" ? error("Can either do manual benders decomposition or use BlockDecomposition.jl!") : nothing
    
    # precalculations
    rf = 1 
    rf_dict = rescaling_factors(rescaling)

    calculate_dependent_values!(network)
    N = nrow(network.buses)
    L = nrow(network.lines)
    T = nrow(network.snapshots) #normally snapshots
    nrow(network.loads_t["p"])!=T ? network.loads_t["p"]=network.loads_t["p_set"] : nothing
    candidates = line_extensions_candidates(network)

    has_reactive_power = !(formulation isa String) && !(formulation <: Union{DCPPowerModel, NFAPowerModel})

    # shortcuts
    buses = network.buses
    generators = network.generators
    links = network.links
    storage_units = network.storage_units
    stores = network.stores
    lines = network.lines

    # indices
    reverse_busidx = rev_idx(buses)
    busidx = idx(buses)
    reverse_lineidx = rev_idx(lines)
    lineidx = idx(lines)

    # iterator bounds
    N_fix_G = sum(((.!generators[:p_nom_extendable]) .& (.!generators[:commitable])))
    N_ext_G = sum(convert(BitArray, generators[:p_nom_extendable]))
    N_com_G = sum(convert(BitArray, generators[:commitable]))
    N_fix_LN = sum((.!lines[:s_nom_extendable]))
    N_ext_LN = sum(.!(.!lines[:s_nom_extendable]))
    N_fix_LK = sum(.!links[:p_nom_extendable])
    N_ext_LK = sum(.!(.!links[:p_nom_extendable]))
    N_fix_SU = sum(.!storage_units[:p_nom_extendable])
    N_ext_SU = sum(.!(.!storage_units[:p_nom_extendable]))
    N_SU = nrow(storage_units)
    N_fix_ST = sum(.!stores[:e_nom_extendable])
    N_ext_ST = sum(.!(.!stores[:e_nom_extendable]))
    N_ST = N_fix_ST + N_ext_ST
    
    if N_fix_LN > 0 && investment_type!="continuous"
        error("Currently all lines have to be extendable for investment type $investment_type!")
    end
    
    if N_com_G > 0
        error("Unit commitment is not implemented yet!")
    end

    sn = snapshot_number
    if sn > 0
        T_params = sn:sn
    else
        T_params = 1:T
    end
    T_params_length = length(T_params)
    
    if blockmodel
        m = BlockModel(solver=solver)
    else
        m = Model(solver=solver)
    end

# --------------------------------------------------------------------------------------------------------
# 2: get base power model from PowerModels
# --------------------------------------------------------------------------------------------------------

    if typeof(formulation) != String
 
        # TODO: unsure, but might be necessary for DC approximation, why does it make such a big difference?
        # if formulation <: GenericPowerModel{T} where T <: pm.AbstractActivePowerFormulation
        #    lines[:r] .= 0.0
        #    lines[:r_pu] .= 0.0
        # end
        
        if formulation <: GenericPowerModel{T} where T <: pm.AbstractBFForm
            post = post_mn_flow_bf
        else
            post = post_mn_flow
        end
        
        network_pm = convert_network(network)
        generic_pm = build_generic_model(network_pm, formulation, post, multinetwork=true, jump_model=m)

    end

# --------------------------------------------------------------------------------------------------------
# 3. add all variables to the Model
# --------------------------------------------------------------------------------------------------------

    println("Add variables")

    if benders != "slave"

        @variable(m, G_p_nom[gr=1:N_ext_G])
        @variable(m, LN_s_nom[l=1:N_ext_LN] >= 0) # for conic constraint inference
        @variable(m, LK_p_nom[l=1:N_ext_LK])
        
        if investment_type == "continuous"
            @variable(m, LN_inv[l=1:N_ext_LN])
        elseif investment_type == "integer"
            @variable(m, LN_inv[l=1:N_ext_LN], Int)
        elseif investment_type == "binary"
            @variable(m, LN_opt[l=1:N_ext_LN], Bin)
            @variable(m, LN_inv[l=1:N_ext_LN])
        elseif investment_type == "integer_bigm"
            @variable(m, LN_opt[l=1:N_ext_LN,c in candidates[l]], Bin)
        else
            error("Investment type $investment_type not defined!")
        end
        
        @variable(m, SU_p_nom[s=1:N_ext_SU])
        @variable(m, ST_e_nom[s=1:N_ext_ST])
    
    end

    if benders != "master"

        @variable(m, G_fix[gr=1:N_fix_G,t=1:T_params_length])
        @variable(m, G_ext[gr=1:N_ext_G,t=1:T_params_length])
        G = [G_fix; G_ext] # G is the concatenated variable array
        
        @variable(m, SU_dispatch_fix[s=1:N_fix_SU,t=1:T_params_length])
        @variable(m, SU_dispatch_ext[s=1:N_ext_SU,t=1:T_params_length])
        SU_dispatch = [SU_dispatch_fix; SU_dispatch_ext]
        
        @variable(m, SU_store_fix[s=1:N_fix_SU,t=1:T_params_length])
        @variable(m, SU_store_ext[s=1:N_ext_SU,t=1:T_params_length])
        SU_store = [SU_store_fix; SU_store_ext]
        
        @variable(m, SU_soc_fix[s=1:N_fix_SU,t=1:T_params_length])
        @variable(m, SU_soc_ext[s=1:N_ext_SU,t=1:T_params_length])
        SU_soc = [SU_soc_fix; SU_soc_ext]
        
        @variable(m, SU_spill_fix[s=1:N_fix_SU,t=1:T_params_length])
        @variable(m, SU_spill_ext[s=1:N_ext_SU,t=1:T_params_length])
        SU_spill = [SU_spill_fix; SU_spill_ext]
        
        @variable(m, ST_dispatch_fix[s=1:N_fix_ST,t=1:T_params_length])
        @variable(m, ST_dispatch_ext[s=1:N_ext_ST,t=1:T_params_length])
        ST_dispatch = [ST_dispatch_fix; ST_dispatch_ext]
        
        @variable(m, ST_store_fix[s=1:N_fix_ST,t=1:T_params_length])
        @variable(m, ST_store_ext[s=1:N_ext_ST,t=1:T_params_length])
        ST_store = [ST_store_fix; ST_store_ext]
        
        @variable(m, ST_soc_fix[s=1:N_fix_ST,t=1:T_params_length])
        @variable(m, ST_soc_ext[s=1:N_ext_ST,t=1:T_params_length])
        ST_soc = [ST_soc_fix; ST_soc_ext]
        
        @variable(m, ST_spill_fix[l=1:N_fix_ST,t=1:T_params_length])
        @variable(m, ST_spill_ext[l=1:N_ext_ST,t=1:T_params_length])
        ST_spill = [ST_spill_fix, ST_spill_ext]
        
        @variable(m, LK_fix[l=1:N_fix_LK,t=1:T_params_length])
        @variable(m, LK_ext[l=1:N_ext_LK,t=1:T_params_length])
        LK = [LK_fix; LK_ext]
        
        if typeof(formulation) != String
            LN_fix = get_LN(network, generic_pm, :p, ext=false)
            LN_ext = get_LN(network, generic_pm, :p, ext=true)
        else
            contains(formulation, "angles") ? @variable(m, THETA[1:N,1:T_params_length]) : nothing
            @variable(m, LN_fix[l=1:N_fix_LN,t=1:T_params_length])
            @variable(m, LN_ext[l=1:N_ext_LN,t=1:T_params_length])
        end
        LN = [LN_fix; LN_ext]

        if has_reactive_power

            # TODO: review if all necessary (potentially merge store/dispatch)
            @variable(m, G_Q_fix[gr=1:N_fix_G,t=1:T_params_length])
            @variable(m, G_Q_ext[gr=1:N_ext_G,t=1:T_params_length])
            G_Q = [G_Q_fix; G_Q_ext]
            
            @variable(m, SU_Q_dispatch_fix[s=1:N_fix_SU,t=1:T_params_length])
            @variable(m, SU_Q_dispatch_ext[s=1:N_ext_SU,t=1:T_params_length])
            SU_Q_dispatch = [SU_Q_dispatch_fix; SU_Q_dispatch_ext]
            
            @variable(m, SU_Q_store_fix[s=1:N_fix_SU,t=1:T_params_length])
            @variable(m, SU_Q_store_ext[s=1:N_ext_SU,t=1:T_params_length])
            SU_Q_store = [SU_Q_store_fix; SU_Q_store_ext]

            @variable(m, ST_Q_dispatch_fix[s=1:N_fix_ST,t=1:T_params_length])
            @variable(m, ST_Q_dispatch_ext[s=1:N_ext_ST,t=1:T_params_length])
            ST_Q_dispatch = [ST_Q_dispatch_fix; ST_Q_dispatch_ext]
            
            @variable(m, ST_Q_store_fix[s=1:N_fix_ST,t=1:T_params_length])
            @variable(m, ST_Q_store_ext[s=1:N_ext_ST,t=1:T_params_length])
            ST_Q_store = [ST_Q_store_fix; ST_Q_store_ext]

            LN_Q_fix = get_LN(network, generic_pm, :q, ext=false)
            LN_Q_ext = get_LN(network, generic_pm, :q, ext=true)
            LN_Q = [LN_Q_fix; LN_Q_ext]

            LN_Q_rev_fix = get_LN(network, generic_pm, :q, ext=false, reverse=true)
            LN_Q_rev_ext = get_LN(network, generic_pm, :q, ext=true, reverse=true)
            LN_Q_rev = [LN_Q_rev_fix; LN_Q_rev_ext]

            LN_rev_fix = get_LN(network, generic_pm, :p, ext=false, reverse=true)
            LN_rev_ext = get_LN(network, generic_pm, :p, ext=true, reverse=true)
            LN_rev = [LN_rev_fix; LN_rev_ext]

        else

            LN_rev = - LN   
        
        end

    end

    if benders == "master"
        @variable(m, ALPHA[g=1:N_cuts]>=0)
    end

    # Determine the order in which the constraints are added to the model and which snapshots belong to it!
    # Useful for plotting constraint matrix.
    # Go through loop only once if full model is built for all snapshots.
    # If subproblems (with subset of snapshots) is built, time indices for variables and parameters are different!

    t_vars_curr = 1
    T_vars_curr = 1

    if !blockstructure && sn==0
        t_vars_curr = T_params_length
        T_vars_curr = T_params
        T_params = [T_params]
    elseif !blockstructure && sn>0
        t_vars_curr = T_params_length
        T_vars_curr = T_params
    end

    for T_params_curr=T_params

        println("Start building model for snapshot(s) $T_params_curr.")

# --------------------------------------------------------------------------------------------------------
# 4. add all generators to the model
# --------------------------------------------------------------------------------------------------------
        
        println("Add generators")

        # set different generator types
        generators = network.generators
        fix_gens_b = ((.!generators[:p_nom_extendable]) .& (.!generators[:commitable]))
        ext_gens_b = convert(BitArray, generators[:p_nom_extendable])
        com_gens_b = convert(BitArray, generators[:commitable])

        p_max_pu = get_switchable_as_dense(network, "generators", "p_max_pu")
        p_min_pu = get_switchable_as_dense(network, "generators", "p_min_pu")

        # add non-extendable generator constraints to the model
        p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=fix_gens_b)
        p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=fix_gens_b)
        filter_timedependent_extremes!(p_max_pu, 0.01)
        filter_timedependent_extremes!(p_min_pu, 0.01)
        p_nom = network.generators[fix_gens_b,:p_nom]

        # slope of power factor
        # TODO: substitute temporary default value
        generators[:pf_max_lag] = 0.8
        generators[:pf_max_lead] = 0.9

        m_lag_ext = pf_slope.(generators[ext_gens_b,:pf_max_lag])
        m_lead_ext = pf_slope.(generators[ext_gens_b,:pf_max_lead])
        m_lag_fix = pf_slope.(generators[fix_gens_b,:pf_max_lag])
        m_lead_fix = pf_slope.(generators[fix_gens_b,:pf_max_lead])

        rf = rf_dict[:bounds_G]

        if benders != "master"
            if blockstructure || sn > 0
                @constraints(m, begin 

                    upper_bounds_G_fix[gr=1:N_fix_G, t=T_params_curr],
                        rf * G_fix[gr,t_vars_curr] <= rf * p_nom[gr] * p_max_pu(t,gr) 
                
                    lower_bounds_G_fix[gr=1:N_fix_G, t=T_params_curr],
                        rf * G_fix[gr,t_vars_curr] >= rf * p_nom[gr] * p_min_pu(t,gr) 
                
                end)
            else
                @constraints(m, begin 
                
                    upper_bounds_G_fix[gr=1:N_fix_G, t=T_params_curr],
                        rf * G_fix[gr,t] <= rf * p_max_pu(t,gr) * p_nom[gr] 
                    
                    lower_bounds_G_fix[gr=1:N_fix_G, t=T_params_curr],
                        rf * G_fix[gr,t] >= rf * p_nom[gr] * p_min_pu(t,gr)  
                
                end)

                has_reactive_power ? @constraints(m, begin

                    upper_bounds_G_Q_fix[gr=1:N_fix_G, t=T_params_curr],
                        rf * G_Q_fix[gr,t] <= rf * p_max_pu(t,gr) * p_nom[gr] * m_lag_fix[gr]  

                    lower_bounds_G_Q_fix[gr=1:N_fix_G, t=T_params_curr],
                        rf * G_Q_fix[gr,t] >= rf * p_nom[gr] * p_max_pu(t,gr) * (-1) * m_lead_fix[gr]
                
                end) : nothing
            end
        end

        # add extendable generator constraints to the model
        p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=ext_gens_b)
        p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=ext_gens_b)
        filter_timedependent_extremes!(p_max_pu, 0.01)
        filter_timedependent_extremes!(p_min_pu, 0.01)
        p_nom_min = network.generators[ext_gens_b,:p_nom_min]
        p_nom_max = network.generators[ext_gens_b,:p_nom_max]

        if benders != "slave" && t_vars_curr==T_params_length
            @constraints(m, begin 
                
                lower_bounds_G_p_nom[gr=1:N_ext_G],
                    G_p_nom[gr] >= p_nom_min[gr]
                
                upper_bounds_G_p_nom[gr=1:N_ext_G],
                    G_p_nom[gr] <= p_nom_max[gr]

            end)
        end

        if benders!="master" && benders!="slave"

            @constraints(m, begin

                upper_bounds_G_ext[t=T_params_curr,gr=1:N_ext_G],
                    rf * G_ext[gr,t] <= rf * p_max_pu(t,gr) * G_p_nom[gr]
 
                lower_bounds_G_ext[t=T_params_curr,gr=1:N_ext_G],
                    rf * G_ext[gr,t] >= rf * p_min_pu(t,gr) * G_p_nom[gr]

            end)
                
            has_reactive_power ? @constraints(m, begin
                
                upper_bounds_G_Q_ext[t=T_params_curr,gr=1:N_ext_G],
                        rf * G_Q_ext[gr,t] <= rf * p_max_pu(t,gr) * G_p_nom[gr] * m_lag_ext[gr]
                    
                lower_bounds_G_Q_ext[t=T_params_curr,gr=1:N_ext_G],
                    rf * G_Q_ext[gr,t] >= rf * p_max_pu(t,gr) * G_p_nom[gr] * (-1) * m_lead_ext[gr]
                                
            end) : nothing

        elseif benders == "slave"

            if blockstructure || sn>0
                @constraints(m, begin
                
                    lower_bounds_G_ext[t=T_params_curr,gr=1:N_ext_G],
                        rf * G_ext[gr,t_vars_curr] >= rf * p_min_pu(t,gr) * generators[ext_gens_b,:p_nom][gr]
                    
                    upper_bounds_G_ext[t=T_params_curr,gr=1:N_ext_G],
                        rf * G_ext[gr,t_vars_curr] <= rf * p_max_pu(t,gr) * generators[ext_gens_b,:p_nom][gr]

                end)
            else
                @constraints(m, begin

                    lower_bounds_G_ext[t=T_params_curr,gr=1:N_ext_G],
                        rf*G_ext[gr,t] >= rf*p_min_pu(t,gr)*generators[ext_gens_b,:p_nom][gr]
                
                    upper_bounds_G_ext[t=T_params_curr,gr=1:N_ext_G],
                        rf*G_ext[gr,t] <= rf*p_max_pu(t,gr)*generators[ext_gens_b,:p_nom][gr]

                end)
            end
        end

        # sort generators the same as corresponding variable array
        generators = [generators[fix_gens_b,:]; generators[ext_gens_b,:] ] 
        fix_gens_b = ((.!generators[:p_nom_extendable]) .& (.!generators[:commitable]))
        ext_gens_b = convert(BitArray, generators[:p_nom_extendable])
        com_gens_b = convert(BitArray, generators[:commitable])

# --------------------------------------------------------------------------------------------------------
# 5. add all lines to the model
# --------------------------------------------------------------------------------------------------------

        println("Add lines")

        # set different lines types
        lines = network.lines
        fix_lines_b = (.!lines[:s_nom_extendable])
        ext_lines_b = .!fix_lines_b

        rf = rf_dict[:bounds_LN]

        # add line constraint for fix lines
        if benders != "master"

            @constraints(m, begin 

                lower_bounds_LN_fix[l=1:N_fix_LN,t=T_params_curr],
                    rf * LN_fix[l,t] >= rf * (-1) * lines[fix_lines_b,:s_max_pu][l] * lines[fix_lines_b,:s_nom][l]  
            
                upper_bounds_LN_fix[l=1:N_fix_LN,t=T_params_curr],
                    rf * LN_fix[l,t] <= rf * lines[fix_lines_b,:s_max_pu][l] * lines[fix_lines_b,:s_nom][l]
          
            end)

        end
        
        if benders != "slave" && t_vars_curr==T_params_length

            @constraints(m, begin        
        
                lower_bounds_LN_s_nom[l=1:N_ext_LN],
                    LN_s_nom[l] >= lines[ext_lines_b,:s_nom_min][l] 
                
                upper_bounds_LN_s_nom[l=1:N_ext_LN],
                    LN_s_nom[l] <= lines[ext_lines_b,:s_nom_max][l]
            
            end)

        end
        
        # add line constraint for extendable lines
        if benders != "master" && benders != "slave"

            @constraints(m, begin

                upper_bounds_LN_ext[t=T_params_curr,l=1:N_ext_LN], 
                    rf * LN_ext[l,t] <=  rf * lines[ext_lines_b,:s_max_pu][l] * LN_s_nom[l]

                lower_bounds_LN_ext[t=T_params_curr,l=1:N_ext_LN], 
                    rf * LN_ext[l,t] >= rf * (-1) * lines[ext_lines_b,:s_max_pu][l] * LN_s_nom[l]

            end)

        elseif benders == "slave"

            if blockstructure || sn>0
                @constraints(m, begin
    
                    upper_bounds_LN_ext[t=T_params_curr,l=1:N_ext_LN], 
                        rf * LN_ext[l,t_vars_curr] <=  
                        rf * lines[:s_max_pu][l] * lines[ext_lines_b,:s_nom][l]
    
                    lower_bounds_LN_ext[t=T_params_curr,l=1:N_ext_LN],
                        rf * LN_ext[l,t_vars_curr] >=
                        rf * (-1) * lines[:s_max_pu][l] * lines[ext_lines_b,:s_nom][l]
    
                end)
            else
                @constraints(m, begin
    
                    upper_bounds_LN_ext[t=T_params_curr,l=1:N_ext_LN], 
                        rf * LN_ext[l,t] <=  
                        rf * lines[:s_max_pu][l] * lines[ext_lines_b,:s_nom][l]
    
                    lower_bounds_LN_ext[t=T_params_curr,l=1:N_ext_LN],
                        rf * LN_ext[l,t] >=
                        rf * (-1) * lines[:s_max_pu][l] * lines[ext_lines_b,:s_nom][l]
    
                end)
            end

        end

        # add logical constraints if applicable
        if benders != "slave"  && t_vars_curr==T_params_length

            if investment_type == "continuous"

                @constraint(m, continuous[l=1:N_ext_LN],
                    LN_s_nom[l] ==
                    ( 1.0 + LN_inv[l] / lines[ext_lines_b,:num_parallel][l] ) * 
                    lines[ext_lines_b,:s_nom][l]
                )
                
            elseif investment_type == "integer"

                @constraint(m, integer[l=1:N_ext_LN],
                    LN_s_nom[l] ==
                    ( 1.0 + LN_inv[l] / lines[ext_lines_b,:num_parallel][l] ) *
                    lines[ext_lines_b,:s_nom][l]
                )
                
            elseif investment_type == "binary"

                bigM_default = 1e4
                bigM = min.(lines[ext_lines_b,:s_nom_max],bigM_default)

                @constraints(m, begin

                    binary1[l=1:N_ext_LN],
                        - bigM[l] * ( 1 - LN_opt[l] ) + lines[ext_lines_b,:s_nom_ext_min][l] <= LN_inv[l]

                    binary2[l=1:N_ext_LN],
                        0 <= LN_inv[l] # no reduction of capacity allowed

                    binary3[l=1:N_ext_LN],
                        bigM[l] * LN_opt[l] >= LN_inv[l]

                    binary4[l=1:N_ext_LN],
                        LN_s_nom[l] == 
                        ( 1.0 + LN_inv[l] / lines[ext_lines_b,:num_parallel][l] ) * lines[ext_lines_b,:s_nom][l]

                end)

            elseif investment_type == "integer_bigm"
                
                @constraint(m, integer_bigm_logic[l=1:N_ext_LN], 
                    sum( LN_opt[l,c] for c in candidates[l] ) == 1.0
                )

                @constraint(m, integer_bigm[l=1:N_ext_LN],
                    LN_s_nom[l] ==
                    ( 1 + sum( c * LN_opt[l,c] for c in candidates[l] ) / lines[ext_lines_b,:num_parallel][l] ) *
                    lines[ext_lines_b,:s_nom][l]
                )
                
            end
        end

        # sort lines the same as corresponding variable array
        lines = [lines[fix_lines_b,:]; lines[ext_lines_b,:]]
        fix_lines_b = (.!lines[:s_nom_extendable])
        ext_lines_b = .!fix_lines_b

# --------------------------------------------------------------------------------------------------------
# 6. add all links to the model
# --------------------------------------------------------------------------------------------------------

    println("Add links")

        rf = rf_dict[:bounds_LK]

        # set different link types
        links = network.links
        fix_links_b = .!links[:p_nom_extendable]
        ext_links_b = .!fix_links_b

        # set constraints for fix links
        if benders != "master"

            if blockstructure || sn > 0
                @constraints(m, begin 

                    lower_bounds_LK_fix[l=1:N_fix_LK,t=T_params_curr],
                        rf * LK_fix[l,t_vars_curr] >= rf * links[fix_links_b, :p_min_pu][l] * links[fix_links_b, :p_nom][l] 
                
                    upper_bounds_LK_fix[l=1:N_fix_LK,t=T_params_curr],
                        rf * LK_fix[l,t_vars_curr] <= rf * links[fix_links_b, :p_max_pu][l] * links[fix_links_b, :p_nom][l]
            
                end)
            else
                @constraints(m, begin 

                    lower_bounds_LK_fix[l=1:N_fix_LK,t=T_params_curr],
                        rf * LK_fix[l,t] >= rf * links[fix_links_b, :p_min_pu][l] * links[fix_links_b, :p_nom][l]
                
                    upper_bounds_LK_fix[l=1:N_fix_LK,t=T_params_curr],
                        rf * LK_fix[l,t] <= rf * links[fix_links_b, :p_max_pu][l] * links[fix_links_b, :p_nom][l]
                
                end)
            end
        end

        if benders != "slave" && t_vars_curr==T_params_length

            @constraints(m, begin 

                lower_bounds_LK_p_nom[l=1:N_ext_LK],
                    LK_p_nom[l] >= links[ext_links_b, :p_nom_min][l]
            
                upper_bounds_LK_p_nom[l=1:N_ext_LK],
                    LK_p_nom[l] <= links[ext_links_b, :p_nom_max][l]
            
            end)
        end

        # set constraints for extendable links
        if benders != "master" && benders != "slave"

            @constraints(m, begin

                lower_bounds_LK_ext[t=T_params_curr,l=1:N_ext_LK],
                    rf * LK_ext[l,t] >= rf * links[ext_links_b, :p_min_pu][l] * LK_p_nom[l]

                upper_bounds_LK_ext[t=T_params_curr,l=1:N_ext_LK],
                    rf * LK_ext[l,t] <= rf * links[ext_links_b, :p_max_pu][l] * LK_p_nom[l]
            end)

        elseif benders == "slave"

            if blockstructure || sn > 0

                @constraint(m, lower_bounds_LK_ext[t=T_params_curr,l=1:N_ext_LK],
                    rf * LK_ext[l,t_vars_curr] >= rf * links[ext_links_b, :p_min_pu][l] * links[ext_links_b,:p_nom][l]
                )

                @constraint(m, upper_bounds_LK_ext[t=T_params_curr,l=1:N_ext_LK],
                    rf * LK_ext[l,t_vars_curr] <= rf * links[ext_links_b, :p_max_pu][l] * links[ext_links_b,:p_nom][l]
                )

            else 

                @constraint(m, lower_bounds_LK_ext[t=T_params_curr,l=1:N_ext_LK],
                    rf * LK_ext[l,t] >= rf * links[ext_links_b, :p_min_pu][l] * links[ext_links_b,:p_nom][l]
                )

                @constraint(m, upper_bounds_LK_ext[t=T_params_curr,l=1:N_ext_LK],
                    rf * LK_ext[l,t] <= rf * links[ext_links_b, :p_max_pu][l] * links[ext_links_b,:p_nom][l]
                )

            end

        end

        # sort links the same as corresponding variable array
        links = [links[fix_links_b,:]; links[ext_links_b,:]]
        fix_links_b = .!links[:p_nom_extendable]
        ext_links_b = .!fix_links_b

# --------------------------------------------------------------------------------------------------------
# 7. define storage_units
# --------------------------------------------------------------------------------------------------------

# WARNING; not implemented for benders and individual snapshots!

        # set different storage_units types
        storage_units = network.storage_units
        fix_sus_b = .!storage_units[:p_nom_extendable]
        ext_sus_b = .!fix_sus_b

        inflow = get_switchable_as_dense(network, "storage_units", "inflow")

        # set constraints for fixed storage_units
        if benders != "master"

            @constraints(m, begin

                lower_bounds_SU_dispatch_fix[s=1:N_fix_SU,t=T_vars_curr], 
                    (0 <=  SU_dispatch_fix[s,t])

                upper_bounds_SU_dispatch_fix[s=1:N_fix_SU,t=T_vars_curr], 
                    (SU_dispatch_fix[s,t] <=
                    (storage_units[fix_sus_b, :p_nom].*storage_units[fix_sus_b, :p_max_pu])[s])
        
                bounds_SU_dispatch_ext[s=1:N_fix_SU,t=T_vars_curr],
                    SU_dispatch_ext[s,t] >= 0
        
                lower_bounds_SU_store_fix[s=1:N_fix_SU,t=T_vars_curr],
                    (0 <=  SU_store_fix[s,t])

                upper_bounds_SU_store_fix[s=1:N_fix_SU,t=T_vars_curr],
                    (SU_store_fix[s,t] <=
                    - (storage_units[fix_sus_b, :p_nom].*storage_units[fix_sus_b, :p_min_pu])[s])
        
                bounds_SU_store_ext[s=1:N_fix_SU,t=T_vars_curr],
                    SU_store_ext[s,t] >= 0
        
                lower_bounds_SU_soc_fix[s=1:N_fix_SU,t=T_vars_curr],
                    0 <= SU_soc_fix[s,t]

                upper_bounds_SU_soc_fix[s=1:N_fix_SU,t=T_vars_curr],
                    SU_soc_fix[s,t] <= 
                    (storage_units[fix_sus_b,:max_hours] .* storage_units[fix_sus_b,:p_nom])[s]
                
                bounds_SU_soc_ext[s=1:N_fix_SU,t=T_vars_curr],
                    SU_soc_ext[s,t] >= 0

            end)
            
            
            if blockstructure || sn>0
                @constraints(m, begin 
                    lower_bounds_SU_spill_fix[s=1:N_fix_SU,t=T_params_curr],
                        0 <=  SU_spill_fix[s,t_vars_curr]
                    
                    lower_bounds_SU_spill_ext[s=1:N_fix_SU,t=T_params_curr],
                        0 <=  SU_spill_ext[s,t_vars_curr]

                    upper_bounds_SU_spill_fix[s=1:N_fix_SU,t=T_params_curr],
                        SU_spill_fix[s,t_vars_curr] <= inflow[:,fix_sus_b][t,s]
                    
                    upper_bounds_SU_spill_ext[s=1:N_fix_SU,t=T_params_curr],
                        SU_spill_ext[s,t_vars_curr] <= inflow[:,ext_sus_b][t,s]
                end)
            else
                @constraints(m, begin 
                    lower_bounds_SU_spill_fix[s=1:N_fix_SU,t=T_params_curr],
                        0 <=  SU_spill_fix[s,t]
                    
                    lower_bounds_SU_spill_ext[s=1:N_fix_SU,t=T_params_curr],
                        0 <=  SU_spill_ext[s,t]

                    upper_bounds_SU_spill_fix[s=1:N_fix_SU,t=T_params_curr],
                        SU_spill_fix[s,t] <= inflow[:,fix_sus_b][t,s]
                    
                    upper_bounds_SU_spill_ext[s=1:N_fix_SU,t=T_params_curr],
                        SU_spill_ext[s,t] <= inflow[:,ext_sus_b][t,s]
                end)
            end
        end
        
        if benders != "slave" && t_vars_curr==T_params_length
            @constraint(m, bounds_SU_p_nom[s=1:N_ext_SU], 
                SU_p_nom[s=1:N_ext_SU] >= 0
            )
        end

        # set constraints for extendable storage_units
        if benders != "master" && benders != "slave"
            @constraints(m, begin
                    su_dispatch_limit[t=T_vars_curr,s=1:N_ext_SU], SU_dispatch_ext[s,t] <= SU_p_nom[s].*storage_units[ext_sus_b, :p_max_pu][s]
                    su_storage_limit[t=T_vars_curr,s=1:N_ext_SU], SU_store_ext[s,t] <= - SU_p_nom[s].*storage_units[ext_sus_b, :p_min_pu][s]
                    su_capacity_limit[t=T_vars_curr,s=1:N_ext_SU], SU_soc_ext[s,t] <= SU_p_nom[s].*storage_units[ext_sus_b, :max_hours][s]
            end)
        elseif benders == "slave"
            @constraints(m, begin
                    su_dispatch_limit[t=T_vars_curr,s=1:N_ext_SU], SU_dispatch_ext[s,t] <= storage_units[ext_sus_b, :p_nom][s].*storage_units[ext_sus_b, :p_max_pu][s]
                    su_storage_limit[t=T_vars_curr,s=1:N_ext_SU], SU_store_ext[s,t] <= - storage_units[ext_sus_b, :p_nom][s].*storage_units[ext_sus_b, :p_min_pu][s]
                    su_capacity_limit[t=T_vars_curr,s=1:N_ext_SU], SU_soc_ext[s,t] <= storage_units[ext_sus_b, :p_nom][s].*storage_units[ext_sus_b, :max_hours][s]
            end)
        end 

        storage_units = [storage_units[fix_sus_b,:]; storage_units[ext_sus_b,:]]
        inflow = [inflow[:,fix_sus_b] inflow[:,ext_sus_b]]

        ext_sus_b = BitArray(storage_units[:p_nom_extendable])

        is_cyclic_i = collect(1:N_SU)[BitArray(storage_units[:cyclic_state_of_charge])]
        not_cyclic_i = collect(1:N_SU)[.!storage_units[:cyclic_state_of_charge]]

        if benders != "master"

            @constraints(m, begin

                    su_logic1[s=is_cyclic_i],
                        SU_soc[s,1] == 
                        (SU_soc[s,T]
                        + storage_units[s,:efficiency_store] * SU_store[s,1]
                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,1]
                        + inflow[1,s] - SU_spill[s,1] )

                    su_logic2[s=not_cyclic_i],
                        SU_soc[s,1] == 
                        (storage_units[s,:state_of_charge_initial]
                        + storage_units[s,:efficiency_store] * SU_store[s,1]
                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,1]
                        + inflow[1,s] - SU_spill[s,1])
        
                    su_logic3[s=1:N_SU,t=2:T],
                        SU_soc[s,t] == 
                        (SU_soc[s,t-1]
                        + storage_units[s,:efficiency_store] * SU_store[s,t]
                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,t]
                        + inflow[t,s] - SU_spill[s,t] )
            end)
        end

# --------------------------------------------------------------------------------------------------------
# 8. define stores
# --------------------------------------------------------------------------------------------------------

        # set different stores types
        stores = network.stores
        fix_stores_b = .!stores[:e_nom_extendable]
        ext_stores_b = .!fix_stores_b

        inflow = get_switchable_as_dense(network, "stores", "inflow")

        #  set constraints for fixed stores
        if benders != "master"

            @constraints(m, begin

                lower_bounds_ST_dispatch_fix[s=1:N_fix_ST,t=T_vars_curr],
                    0 <=  ST_dispatch_fix[s,t]

                upper_bounds_ST_dispatch_fix[s=1:N_fix_ST,t=T_vars_curr],
                    (ST_dispatch_fix[s,t] <=
                    (stores[fix_stores_b, :e_nom].*stores[fix_stores_b, :e_max_pu])[s])

                bounds_ST_dispatch_ext[s=1:N_fix_ST,t=T_vars_curr],
                    ST_dispatch_ext[s,t] >= 0

                lower_bounds_ST_store_fix[s=1:N_fix_ST,t=T_vars_curr],
                    0 <=  ST_store_fix[s,t]

                upper_bounds_ST_store_fix[s=1:N_fix_ST,t=T_vars_curr],
                    (ST_store_fix[s,t] <=
                    - (stores[fix_stores_b, :e_nom].*stores[fix_stores_b, :e_min_pu])[s])

                bounds_ST_store_ext[s=1:N_fix_ST,t=T_vars_curr],
                    ST_store_ext[s,t] >= 0

                lower_bounds_ST_soc_fix[s=1:N_fix_ST,t=T_vars_curr],
                    0 <= ST_soc_fix[s,t]

                upper_bounds_ST_soc_fix[s=1:N_fix_ST,t=T_vars_curr],
                    ST_soc_fix[s,t] <= 
                    (stores[fix_stores_b,:max_hours] .* stores[fix_stores_b,:e_nom])[s]

                bounds_ST_soc_ext[s=1:N_fix_ST,t=T_vars_curr],
                    ST_soc_ext[s,t] >= 0

            end)
            
            
            if blockstructure || sn>0
                @constraints(m, begin 
                    lower_bounds_ST_spill_fix[s=1:N_fix_ST,t=T_params_curr],
                        0 <=  ST_spill_fix[s,t_vars_curr]
    
                    lower_bounds_ST_spill_ext[s=1:N_fix_ST,t=T_params_curr],
                        0 <=  ST_spill_ext[s,t_vars_curr]

                    upper_bounds_ST_spill_fix[s=1:N_fix_ST,t=T_params_curr],
                        ST_spill_fix[s,t_vars_curr] <= 
                        inflow[:,fix_stores_b][s,t]
    
                    upper_bounds_ST_spill_ext[s=1:N_fix_ST,t=T_params_curr],
                        ST_spill_ext[s,t_vars_curr] <= 
                        inflow[:,ext_stores_b][s,t]
                end)
            else 
                @constraints(m, begin 
                    lower_bounds_ST_spill_fix[s=1:N_fix_ST,t=T_params_curr],
                        0 <=  ST_spill_fix[s,t]
    
                    lower_bounds_ST_spill_ext[s=1:N_fix_ST,t=T_params_curr],
                        0 <=  ST_spill_ext[s,t]
                  
                    upper_bounds_ST_spill_fix[s=1:N_fix_ST,t=T_params_curr],
                        ST_spill_fix[s,t] <= 
                        inflow[:,fix_stores_b][s,t]
    
                    upper_bounds_ST_spill_ext[s=1:N_fix_ST,t=T_params_curr],
                        ST_spill_ext[s,t] <= 
                        inflow[:,ext_stores_b][s,t]
                end)
            end
        end
        
        if benders != "slave" && t_vars_curr==T_params_length
            @constraint(m, bounds_ST_e_nom[s=1:N_ext_ST],
                ST_e_nom[s] >= 0
            )
        end

        # set constraints for extendable stores
        if benders != "master" && benders != "slave"
            @constraints(m, begin
                    st_dispatch_limit[t=T_vars_curr,s=1:N_ext_ST], 
                        ST_dispatch_ext[s,t] <= ST_e_nom[s].*stores[ext_stores_b, :e_max_pu][s]
                    st_storage_limit[t=T_vars_curr,s=1:N_ext_ST], 
                        ST_store_ext[s,t] <= - ST_e_nom[s].*stores[ext_stores_b, :e_min_pu][s]
                    st_capacity_limit[t=T_vars_curr,s=1:N_ext_ST], 
                        ST_soc_ext[s,t] <= ST_e_nom[s].*stores[ext_stores_b, :max_hours][s]
            end)
        elseif benders == "slave"
            @constraints(m, begin
                    st_dispatch_limit[t=T_vars_curr,s=1:N_ext_ST],
                        ST_dispatch_ext[s,t] <= 
                        stores[ext_stores_b, :e_nom][s].*stores[ext_stores_b, :e_max_pu][s]
                    st_storage_limit[t=T_vars_curr,s=1:N_ext_ST],
                        ST_store_ext[s,t] <= 
                        - stores[ext_stores_b, :e_nom][s].*stores[ext_stores_b, :e_min_pu][s]
                    st_capacity_limit[t=T_vars_curr,s=1:N_ext_ST],
                        ST_soc_ext[s,t] <= 
                        stores[ext_stores_b, :e_nom][s].*stores[ext_stores_b, :max_hours][s]
            end)
        end

        # set charging constraints
        stores = [stores[fix_stores_b,:]; stores[ext_stores_b,:]]
        inflow = [inflow[:,fix_stores_b] inflow[:,ext_stores_b]]

        ext_stores_b = BitArray(stores[:e_nom_extendable])

        is_cyclic_i = collect(1:N_ST)[BitArray(stores[:cyclic_state_of_charge])]
        not_cyclic_i = collect(1:N_ST)[.!stores[:cyclic_state_of_charge]]

        if benders != "master"

            @constraints(m, begin

                    st_logic1[s=is_cyclic_i,t=1], 
                        ST_soc[s,t] == 
                        (ST_soc[s,T]
                        + stores[s,:efficiency_store] * ST_store[s,t]
                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]
                        + inflow[t,s] - ST_spill[s,t])
        
                    st_logic2[s=not_cyclic_i,t=1], 
                        ST_soc[s,t] == 
                        (stores[s,:state_of_charge_initial]
                        + stores[s,:efficiency_store] * ST_store[s,t]
                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]
                        + inflow[t,s] - ST_spill[s,t])
        
                    st_logic3[s=is_cyclic_i,t=2:T], 
                        ST_soc[s,t] == 
                        (ST_soc[s,t-1]
                        + stores[s,:efficiency_store] * ST_store[s,t]
                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]
                        + inflow[t,s] - ST_spill[s,t])

            end)
        end

# --------------------------------------------------------------------------------------------------------
# 9. add shunt impedances
# --------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------
# 10. power flow formulations
# --------------------------------------------------------------------------------------------------------

        println("Add flow constraints")

        rf = rf_dict[:flows]

        if benders != "master" 
            
            if typeof(formulation) != String

                # load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                # TODO: add shunts
              
                @constraint(m, nodal_p[t=T_vars_curr,n=1:N],
    
                        sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                        + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                        .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                        + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                        - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                        - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                        + sum(LN_rev[findin(lines[:bus1], [reverse_busidx[n]]),t]) 
                )


                if has_reactive_power
                    @constraint(m, nodal_q[t=T_vars_curr,n=1:N],
        
                            sum(G_Q[findin(generators[:bus], [reverse_busidx[n]]), t])
                            + sum(SU_Q_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                            #- row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                            - sum(SU_Q_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                            == sum(LN_Q[findin(lines[:bus0], [reverse_busidx[n]]),t])
                            + sum(LN_Q_rev[findin(lines[:bus1], [reverse_busidx[n]]),t]) 
                    )
                end

                # add cm upper bounds dependant on line capacity

                if typeof(formulation) != String && formulation <: Union{SOCBFPowerModel, SOCBFConicPowerModel, QCWRPowerModel, QCWRTriPowerModel}

                    acc = pm.ref(generic_pm, 1, :bus)

                    # TODO: add option in function parameters
                    convex_relaxation = true

                    ub_ext = Array{Any}(N_ext_LN)

                    if convex_relaxation
                        for l=1:N_ext_LN

                            point_at_lim(x) = ( x, (lines[ext_lines_b,:s_max_pu][l] * x / acc[busidx[lines[l,:bus0]]]["vmin"])^2) 
                            # TODO: access vminpu from PSA once default set

                            min_point = point_at_lim(lines[ext_lines_b,:s_nom_min][l])
                            max_point = point_at_lim(lines[ext_lines_b,:s_nom_max][l])

                            slope = get_slope(min_point, max_point)
                            intercept = get_intercept(min_point, max_point)

                            ub_ext[l] = slope * LN_s_nom[l] + intercept

                        end
                    else
                        for l=1:N_ext_LN
                            ub_ext[l] = (lines[ext_lines_b,:s_max_pu][l] * LN_s_nom[l] / acc[busidx[lines[l,:bus0]]]["vmin"])^2
                        end
                    end

                end

                if typeof(formulation) != String && formulation <: Union{QCWRPowerModel, QCWRTriPowerModel}
                    
                    @constraint(generic_pm.model,
                        ub_cm_ext[l=1:N_ext_LN,t=T_vars_curr],
                        pm.var(generic_pm, t, 1)[:cm][(busidx[lines[l,:bus0]],busidx[lines[l,:bus1]])] <= ub_ext[l]
                    ) 

                    @constraint(generic_pm.model,
                        ub_cm_fix[l=1:N_fix_LN,t=T_vars_curr],
                        pm.var(generic_pm, t, 1)[:cm][(busidx[lines[l,:bus0]],busidx[lines[l,:bus1]])] <= 
                        (lines[fix_lines_b,:s_max_pu][l] * lines[fix_lines_b,:s_nom][l] / acc[busidx[lines[l,:bus0]]]["vmin"])^2
                    ) # TODO: access vminpu from PSA once default set

                elseif formulation <: Union{SOCBFPowerModel, SOCBFConicPowerModel}

                    @constraint(generic_pm.model,
                        ub_cm_ext[l=1:N_ext_LN,t=T_vars_curr],
                        pm.var(generic_pm, t, 1,:cm)[lineidx[lines[l,:name]]] <= ub_ext[l]
                    )
                    
                    @constraint(generic_pm.model,
                        ub_cm_fix[l=1:N_fix_LN,t=T_vars_curr],
                        pm.var(generic_pm, t, 1,:cm)[lineidx[lines[l,:name]]] <= 
                        (lines[fix_lines_b,:s_max_pu][l] * lines[fix_lines_b,:s_nom][l] / acc[busidx[lines[l,:bus0]]]["vmin"])^2
                    ) # TODO: access vminpu from PSA once default set

                end

                if typeof(formulation) != String && formulation <: Union{SOCWRPowerModel, SOCWRConicPowerModel, SOCBFPowerModel, SOCBFConicPowerModel, SDPWRMPowerModel, SparseSDPWRMPowerModel, QCWRPowerModel, QCWRTriPowerModel}
                    
                    ids_ext = map(x -> reverse_lineidx[x], lines[ext_lines_b,:][:name])
                    ids_fix = map(x -> reverse_lineidx[x], lines[fix_lines_b,:][:name])
                    
                    for (n, netw) in nws(generic_pm)

                        l=1
                        for i in ids_fix
                            bound = lines[fix_lines_b,:s_max_pu][l] * lines[fix_lines_b,:s_nom][l]
                            constraint_thermal_limit_from_fix(generic_pm, i, bound, nw=n)
                            constraint_thermal_limit_to_fix(generic_pm, i, bound, nw=n)
                            l+=1
                        end

                        l = 1
                        for i in ids_ext
                            bound = lines[ext_lines_b,:s_max_pu][l]*LN_s_nom[l]
                            constraint_thermal_limit_from_ext(generic_pm, i, bound, nw=n)
                            constraint_thermal_limit_to_ext(generic_pm, i, bound, nw=n)
                            l+=1
                        end

                    end
                end
                

            elseif formulation == "angles_linear"

                # load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                if sn > 0

                    @constraint(m, nodal[t=T_vars_curr,n=1:N],
    
                            sum(G[findin(generators[:bus], [reverse_busidx[n]]), t_vars_curr])
                            + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t_vars_curr])
                            + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t_vars_curr])
    
                            - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                            - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t_vars_curr])
                            - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t_vars_curr])
    
                            == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t_vars_curr])
                            - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t_vars_curr]) 
                    )

                else

                    @constraint(m, nodal[t=T_vars_curr,n=1:N],
    
                            sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                            + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                            + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
    
                            - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                            - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                            - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
    
                            == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                            - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) 
                    )
                end

                sn > 0 ? cnt = t_vars_curr : cnt = T_vars_curr

                @constraint(m, flows[t=cnt,l=1:L], 
                    rf * LN[l, t] == 
                    rf * lines[:x_pu][l]^(-1) *     
                    ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                )

                @constraint(m, slack[t=cnt], THETA[1,t] == 0 )

            elseif formulation == "angles_linear_integer_bigm"

                # get bigm parameters for each line
                bigm_upper = bigm(:flows_upper, network)
                bigm_lower = bigm(:flows_lower, network)

                # load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                if benders!="master" && benders!="slave"

                    @constraint(m, nodal[t=T_vars_curr,n=1:N], (

                        sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                        + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                        + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                        - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                        - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                        - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) ))

                    @constraint(m, flows_upper[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=T_vars_curr],
                        rf * (
                            ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) * 
                            ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                            - LN[l,t]
                        ) >= rf * ( LN_opt[l,c] - 1 ) * bigm_upper[l]
                    )
                    
                    @constraint(m, flows_lower[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=T_vars_curr],
                        rf * (
                            ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) * 
                            ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                            - LN[l,t]
                        ) <= rf * ( 1 - LN_opt[l,c] ) * bigm_lower[l]    
                    )

                    @constraint(m, flows_fix[l=1:(sum(fix_lines_b)), t=T_vars_curr],

                        rf * LN[l, t] == 
                        rf * lines[:x_pu][l]^(-1) * ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) 
                    
                    )   

                    @constraint(m, slack[t=T_vars_curr], THETA[1,t] == 0)

                elseif benders=="slave"

                    if blockstructure || sn>0

                        @constraint(m, nodal[t=T_vars_curr,n=1:N], (

                            sum(G[findin(generators[:bus], [reverse_busidx[n]]), t_vars_curr])
                            + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                                .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t_vars_curr])
                            + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t_vars_curr])

                            - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                            - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t_vars_curr])
                            - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t_vars_curr])

                            == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t_vars_curr])
                            - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t_vars_curr]) ))

                        @constraint(m, flows_upper[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=T_vars_curr],
                            rf * (
                                ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) * 
                                ( THETA[busidx[lines[:bus0][l]], t_vars_curr] - THETA[busidx[lines[:bus1][l]], t_vars_curr] ) 
                                - LN[l,t_vars_curr]
                            ) >= rf * bigm_upper[l] * (c == 0 ? 0 : -1)
                        )
        
                        @constraint(m, flows_lower[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=T_vars_curr],
                            rf * (
                                ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) * 
                                ( THETA[busidx[lines[:bus0][l]], t_vars_curr] - THETA[busidx[lines[:bus1][l]], t_vars_curr] ) 
                                - LN[l,t_vars_curr]
                            ) <= rf * bigm_lower[l] * (c == 0 ? 0 : 1)
                        )

                        @constraint(m, flows_fix[l=1:(sum(fix_lines_b)), t=T_vars_curr],
                            rf * LN[l, t_vars_curr] == 
                            rf * lines[:x_pu][l]^(-1) * ( THETA[busidx[lines[:bus0][l]], t_vars_curr] - THETA[busidx[lines[:bus1][l]], t_vars_curr] ) 
                        ) 
                        
                        @constraint(m, slack[t=T_vars_curr], THETA[1,t_vars_curr] == 0)

                    else

                        @constraint(m, nodal[t=T_vars_curr,n=1:N], (

                            sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                            + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                                .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                            + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                            - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                            - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                            - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                            == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                            - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) ))

                        @constraint(m, flows_upper[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=T_params_curr],
                            rf * (
                                ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) * 
                                ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                                - LN[l,t]
                            ) >= rf * bigm_upper[l] * (c == 0 ? 0 : -1)
                        )
        
                        @constraint(m, flows_lower[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=T_params_curr],
                            rf * (
                                ( 1 + c / lines[:num_parallel][l] ) *lines[:x_pu][l]^(-1) *
                                ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                                - LN[l,t]
                            ) <= rf * bigm_lower[l] * (c == 0 ? 0 : 1)
                        )

                        @constraint(m, flows_nonext[l=1:(sum(fix_lines_b)), t=T_params_curr],
                            LN[l, t] == lines[:x_pu][l]^(-1) *
                            ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) 
                        )  

                        @constraint(m, slack[t=T_params_curr], THETA[1,t] == 0)

                    end
                end   

            elseif formulation == "angles_bilinear"

                # cannot be solved by Gurobi, needs Ipopt!
                # requires starting reactance and line capacity!

                # load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                @constraint(m, nodal[t=T_vars_curr,n=1:N], (

                        sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                        + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                                .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                        + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                        - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                        - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                        - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) )        
                )

                @NLconstraint(m, flows_ext[t=T_vars_curr,l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b))],
                    LN[l, t] ==  
                    (1 + LN_inv[l] / lines[:num_parallel][l]) * lines[:x_pu][l]^(-1) *
                    (THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) 
                )

                @constraint(m, flows_nonext[t=T_vars_curr,l=1:(sum(fix_lines_b))],
                    LN[l, t] == lines[:x_pu][l]^(-1) *
                    ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) 
                )                                                   

                @constraint(m, slack[t=T_vars_curr], THETA[1,t] == 0)

            elseif formulation == "kirchhoff_linear"

                #load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                @constraint(m, balance[t=T_vars_curr,n=1:N], (
                    sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                    + sum(LN[ findin(lines[:bus1], [reverse_busidx[n]]) ,t])
                    + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                        .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                    + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                    - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                    - sum(LN[ findin(lines[:bus0], [reverse_busidx[n]]) ,t])
                    - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                    - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                    == 0 )
                )

                # Might be nessecary to loop over all subgraphs as
                # for (sn, sub) in enumerate(weakly_connected_components(g))
                #     g_sub = induced_subgraph(g, sub)[1]

                append_idx_col!(lines)
                (branches, var, attribute) = (lines, LN, :x)
                cycles = get_cycles(network)
                if ndims(cycles)<2
                    cycles = [cycle for cycle in cycles if length(cycle)>2]
                else
                    cycles = [cycles[i,:] for i in 1:size(cycles)[1]]
                end
                if length(cycles)>0
                    cycles_branch = Array{Int64,1}[]
                    directions = Array{Float64,1}[]
                    for cyc=1:length(cycles)
                        push!(cycles_branch,Int64[])
                        push!(directions,Float64[])
                        for bus=1:length(cycles[cyc])
                            bus0 = cycles[cyc][bus]
                            bus1 = cycles[cyc][(bus)%length(cycles[cyc])+1]
                            try
                                push!(cycles_branch[cyc],branches[((branches[:bus0].==reverse_busidx[bus0])
                                            .&(branches[:bus1].==reverse_busidx[bus1])),:idx][1] )
                                push!(directions[cyc], 1.)
                            catch y
                                if isa(y, BoundsError)
                                    push!(cycles_branch[cyc], branches[((branches[:bus0].==reverse_busidx[bus1])
                                                    .&(branches[:bus1].==reverse_busidx[bus0])),:idx][1] )
                                    push!(directions[cyc], -1.)
                                else
                                    return y
                                end
                            end
                        end
                    end
                    if attribute==:x
                        @constraint(m, line_cycle_constraint[t=T_vars_curr,c=1:length(cycles_branch)] ,
                                dot(directions[c] .* lines[cycles_branch[c], :x_pu],
                                    LN[cycles_branch[c],t]) == 0
                        )
                    # elseif attribute==:r
                    #     @constraint(m, link_cycle_constraint[c=1:length(cycles_branch), t=T_vars_curr] ,
                    #             dot(directions[c] .* links[cycles_branch[c], :r]/380. , LK[cycles_branch[c],t]) == 0)
                    end
                end

            elseif formulation == "kirchhoff_bilinear"

                # cannot be solved by Gurobi, needs Ipopt!
                # requires starting reactance and line capacity!

                #load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                @constraint(m, balance[t=T_vars_curr,n=1:N], (

                    sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                    + sum(LN[ findin(lines[:bus1], [reverse_busidx[n]]) ,t])
                    + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                        .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                    + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                    - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                    - sum(LN[ findin(lines[:bus0], [reverse_busidx[n]]) ,t])
                    - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                    - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                    == 0 )
                )

                # Might be nessecary to loop over all subgraphs as
                # for (sn, sub) in enumerate(weakly_connected_components(g))
                #     g_sub = induced_subgraph(g, sub)[1]

                append_idx_col!(lines)
                (branches, var, attribute) = (lines, LN, :x)
                cycles = get_cycles(network)
                if ndims(cycles)<2
                    cycles = [cycle for cycle in cycles if length(cycle)>2]
                else
                    cycles = [cycles[i,:] for i in 1:size(cycles)[1]]
                end
                if length(cycles)>0
                    cycles_branch = Array{Int64,1}[]
                    directions = Array{Float64,1}[]
                    for cyc=1:length(cycles)
                        push!(cycles_branch,Int64[])
                        push!(directions,Float64[])
                        for bus=1:length(cycles[cyc])
                            bus0 = cycles[cyc][bus]
                            bus1 = cycles[cyc][(bus)%length(cycles[cyc])+1]
                            try
                                push!(cycles_branch[cyc],branches[((branches[:bus0].==reverse_busidx[bus0])
                                            .&(branches[:bus1].==reverse_busidx[bus1])),:idx][1] )
                                push!(directions[cyc], 1.)
                            catch y
                                if isa(y, BoundsError)
                                    push!(cycles_branch[cyc], branches[((branches[:bus0].==reverse_busidx[bus1])
                                                    .&(branches[:bus1].==reverse_busidx[bus0])),:idx][1] )
                                    push!(directions[cyc], -1.)
                                else
                                    return y
                                end
                            end
                        end
                    end
                    if attribute==:x
                        @NLconstraint(m, line_cycle_constraint[t=T_vars_curr,c=1:length(cycles_branch)] ,
                                sum(      directions[c][l] 
                                        * lines[cycles_branch[c], :x_pu][l] 
                                        * (1+LN_inv[cycles_branch[c]][l] / lines[cycles_branch[c], :num_parallel][l])^(-1)
                                        #--* (1+LN_inv[cycles_branch[c]][l])^(-1)
                                        * LN[cycles_branch[c],t][l]
                                for l=1:length(directions[c])
                                )               
                                    == 0)
                    # elseif attribute==:r
                    #     @constraint(m, link_cycle_constraint[c=1:length(cycles_branch), t=T_vars_curr] ,
                    #             dot(directions[c] .* links[cycles_branch[c], :r]/380. , LK[cycles_branch[c],t]) == 0)
                    end
                end

            elseif formulation == "ptdf"

                ptdf = ptdf_matrix(network)

                #load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                @constraint(m, flows[t=T_vars_curr,l=1:L],
                    sum( ptdf[l,n]
                    * (   sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                        + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                        + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                        - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                        - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
                        )
                    for n in 1:N
                    ) == LN[l,t] 
                )

                @constraint(m, balance[t=T_vars_curr],
                    sum( (   sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                        + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                        + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                        - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                        - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
                        ) for n in 1:N
                    ) == 0 
                )

            else
                error("The formulation $formulation is not implemented.")
            end
        end
    

# --------------------------------------------------------------------------------------------------------
# 11. set global_constraints
# --------------------------------------------------------------------------------------------------------

        println("Add global constraints")

        # carbon constraint
        ## carbon emissions <= emission limit
        if t_vars_curr==T_params_length

            carrier_index(carrier) = findin(generators[:carrier], carrier)
    
            if benders != "master" && sn==0
                if nrow(network.global_constraints)>0 && in("co2_limit", network.global_constraints[:name])
                    co2_limit = network.global_constraints[network.global_constraints[:name].=="co2_limit", :constant]
                    println("CO2_limit is $(co2_limit[1]) t")
                    nonnull_carriers = network.carriers[network.carriers[:co2_emissions].!=0, :][:name]
                    @constraint(m, co2limit, sum(sum(network.snapshots[:weightings][t]*dot(1./generators[carrier_index(nonnull_carriers) , :efficiency],
                                G[carrier_index(nonnull_carriers),t]) for t=1:T)
                                * select_names(network.carriers, [carrier])[:co2_emissions]
                                for carrier in network.carriers[:name]) .<=  co2_limit)
                end
            end
    
            # limit on transmission expansion volume
            ## sum of capacity expansion times length over all lines <= limit in MWkm unit
            if benders != "slave"
                if nrow(network.global_constraints)>0 && in(true, network.lines[:s_nom_extendable]) && in("mwkm_limit", network.global_constraints[:name])
                    mwkm_limit = network.global_constraints[network.global_constraints[:name].=="mwkm_limit", :constant]
                    #println("Line expansion limit is $(mwkm_limit[1]) times current MWkm")
                    @constraint(m, mwkmlimit, 
                        dot(LN_s_nom,lines[:length]) <= mwkm_limit[1] * dot(lines[ext_lines_b,:s_nom],lines[:length])
                    )
                end
            end
    
            # renewable energy target
            ## sum of renewable generation =/>= percentage * sum of total load
            if benders != "master" && sn==0
                if nrow(network.global_constraints)>0 && in("restarget", network.global_constraints[:name])
                    restarget = network.global_constraints[network.global_constraints[:name].=="restarget", :constant]
                    #println("Target share of renewable energy is $(restarget[1]*100) %")
                    null_carriers = network.carriers[network.carriers[:co2_emissions].==0,:][:name]
                    @constraint(m, restarget,
                        sum(sum(network.snapshots[:weightings][t]*G[carrier_index(null_carriers),t] for t=1:T))
                        .>= restarget * sum(network.snapshots[:weightings][t]*sum(convert(Array,network.loads_t["p_set"][t,2:end])) for t=1:T)
                    )
                end
            end
    
            # approximate renewable energy target
            ## approximation based on investment variables only!
            ## idea: curtailment + (variable) renewable generation >= percentage * sum of total load
            ## prerequisite: curtailment is low
            if benders != "slave"
                if nrow(network.global_constraints)>0 && in("approx_restarget", network.global_constraints[:name])
    
                    N_loads = size(network.loads_t["p_set"])[2]
                    approx_restarget = network.global_constraints[network.global_constraints[:name].=="approx_restarget", :constant]
    
                    # needed to take out biomass since it is not time dependent
                    null_carriers = network.carriers[(network.carriers[:co2_emissions].==0) .& (network.carriers[:name] .!= "biomass"),:][:name]
                    
                    ren_gens_b = [in(i,carrier_index(null_carriers)) ? true : false for i=1:size(generators)[1]]
                    fix_ren_gens_b = .!generators[:p_nom_extendable][ren_gens_b]
                    ext_ren_gens_b = .!fix_ren_gens_b
                    
                    ren_gens_b_orig = [in(i,findin(network.generators[:carrier], null_carriers)) ? true : false for i=1:size(network.generators)[1]]
                    fix_ren_gens_b_orig = .!network.generators[:p_nom_extendable][ren_gens_b_orig]
                    ext_ren_gens_b_orig = .!fix_ren_gens_b_orig
                    
                    fix_ren_gens = generators[fix_gens_b .& ren_gens_b,:]
                    ext_ren_gens = generators[ext_gens_b .& ren_gens_b,:]
                    
                    def_p_max_pu_ext = 8760.0*ext_ren_gens[:p_max_pu]
                    def_p_max_pu_fix = 8760.0*fix_ren_gens[:p_max_pu]
                    
                    exist_fix_ren_gens = maximum(fix_ren_gens_b_orig)

                    p_max_pu_full = network.generators_t["p_max_pu"][:,2:end]
                    exist_fix_ren_gens ? p_max_pu_fix = convert(Array,p_max_pu_full[:,fix_ren_gens_b_orig]) : nothing
                    p_max_pu_ext = convert(Array,p_max_pu_full[:,ext_ren_gens_b_orig])
    
                    loc_fix = findin(fix_ren_gens[:name], string.(p_max_pu_full.colindex.names))
                    loc_ext = findin(ext_ren_gens[:name], string.(p_max_pu_full.colindex.names))
    
                    loc_fix_b = [in(i,loc_fix) ? true : false for i=1:length(def_p_max_pu_fix)]
                    loc_ext_b = [in(i,loc_ext) ? true : false for i=1:length(def_p_max_pu_ext)]
                    
                    exist_fix_ren_gens ? sum_of_p_max_pu_fix = sum(network.snapshots[:weightings][t]*p_max_pu_fix[t,:] for t=1:T)  : nothing
                    sum_of_p_max_pu_ext = sum(network.snapshots[:weightings][t]*p_max_pu_ext[t,:] for t=1:T)
    
                    exist_fix_ren_gens ? def_p_max_pu_fix[loc_fix_b] .= sum_of_p_max_pu_fix : nothing
                    def_p_max_pu_ext[loc_ext_b] .= sum_of_p_max_pu_ext

                    rf = rf_dict[:approx_restarget]

                    if exist_fix_ren_gens

                        @constraint(m, approx_restarget,
                            rf * (dot(def_p_max_pu_fix,fix_ren_gens[:p_nom]) + dot(def_p_max_pu_ext,G_p_nom))
                            >= rf * approx_restarget[1] * sum(network.snapshots[:weightings][t]*network.loads_t["p_set"][t,n] for t=1:T for n=2:N_loads)
                        )
                        
                    else
                        @constraint(m, approx_restarget,
                            rf * dot(def_p_max_pu_ext,G_p_nom)
                            >= rf * approx_restarget[1] * sum(network.snapshots[:weightings][t]*network.loads_t["p_set"][t,n] for t=1:T for n=2:N_loads)
                        )
                    end
                end
            end
        end
        

# --------------------------------------------------------------------------------------------------------
# 12. set objective function
# --------------------------------------------------------------------------------------------------------
    
        println("Add objective")

        if t_vars_curr==T_params_length

            rf = rf_dict[:objective]

            if benders!="master" && benders!="slave"
    
                @objective(m, Min,
                    rf * (
                        sum(network.snapshots[:weightings][t] * dot(generators[:marginal_cost], G[:,t]) for t=1:T_params_length)
                        + dot(generators[ext_gens_b,:capital_cost], G_p_nom[:] )
                        + dot(generators[fix_gens_b,:capital_cost], generators[fix_gens_b,:p_nom])
    
                        + dot(lines[ext_lines_b,:capital_cost], LN_s_nom[:])
                        + dot(lines[fix_lines_b,:capital_cost], lines[fix_lines_b,:s_nom])
    
                        + dot(links[ext_links_b,:capital_cost], LK_p_nom[:])
                        + dot(links[fix_links_b,:capital_cost], links[fix_links_b,:p_nom])
    
                        + sum(network.snapshots[:weightings][t] * dot(storage_units[:marginal_cost], SU_dispatch[:,t]) for t=1:T_params_length)
                        + dot(storage_units[ext_sus_b, :capital_cost], SU_p_nom[:])
                        + dot(storage_units[fix_sus_b,:capital_cost], storage_units[fix_sus_b,:p_nom])
    
                        + sum(network.snapshots[:weightings][t] * dot(stores[:marginal_cost], ST_dispatch[:,t]) for t=1:T_params_length)
                        + dot(stores[ext_stores_b, :capital_cost], ST_e_nom[:])
                        + dot(stores[fix_stores_b,:capital_cost], stores[fix_stores_b,:e_nom])
                    )
                )
    
            elseif benders == "master"
    
                @objective(m, Min,
                    rf * (
                    dot(generators[ext_gens_b,:capital_cost], G_p_nom[:] )
                    + dot(generators[fix_gens_b,:capital_cost], generators[fix_gens_b,:p_nom])
            
                    + dot(lines[ext_lines_b,:capital_cost], LN_s_nom[:])
                    + dot(lines[fix_lines_b,:capital_cost], lines[fix_lines_b,:s_nom])
            
                    + dot(links[ext_links_b,:capital_cost], LK_p_nom[:])
                    + dot(links[fix_links_b,:capital_cost], links[fix_links_b,:p_nom])
            
                    + dot(storage_units[ext_sus_b, :capital_cost], SU_p_nom[:])
                    + dot(storage_units[fix_sus_b,:capital_cost], storage_units[fix_sus_b,:p_nom])
            
                    + dot(stores[ext_stores_b, :capital_cost], ST_e_nom[:])
                    + dot(stores[fix_stores_b,:capital_cost], stores[fix_stores_b,:e_nom])
                    )
                    + sum( ALPHA[g] for g=1:N_cuts )
                )

            elseif benders == "slave"
                
                if sn>0
                    @objective(m, Min,
                        rf * (
                              sum( network.snapshots[:weightings][t] * dot(generators[:marginal_cost], G[:,t_vars_curr]) for t=T_vars_curr )
                            + sum( network.snapshots[:weightings][t] * dot(storage_units[:marginal_cost], SU_dispatch[:,t_vars_curr]) for t=T_vars_curr )
                            + sum( network.snapshots[:weightings][t] * dot(stores[:marginal_cost], ST_dispatch[:,t_vars_curr]) for t=T_vars_curr )
                        )
                    )
                else
                    @objective(m, Min,
                        rf * (
                              sum( network.snapshots[:weightings][t] * dot(generators[:marginal_cost], G[:,t]) for t=T_vars_curr )
                            + sum( network.snapshots[:weightings][t] * dot(storage_units[:marginal_cost], SU_dispatch[:,t]) for t=T_vars_curr )
                            + sum( network.snapshots[:weightings][t] * dot(stores[:marginal_cost], ST_dispatch[:,t]) for t=T_vars_curr )
                        )
                    )
                end
    
            else
                error()
            end
        end

        t_vars_curr += 1
        T_vars_curr = t_vars_curr

    end # for loop for problem building

    return m

end
