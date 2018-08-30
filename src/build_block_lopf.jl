using BlockDecomposition

# TODO:
# 1. what decision variables go in master, which in sub? what are the problem sizes?
# --- option A: all investment variables in master
# --- option B: only line investment in master
# 2. add decision variable for reactance in master problem (requires exisitng line)
# --- LN_x_pu[l] = ( network.lines[:s_nom][l]  * network.lines[:x_pu][l] ) / LN_s_nom[l]
# --- added to the master problem this is non-linear! - But linear with susceptance!
# --- what would happen if x was updated everytime a master solution is handed to the master problem?
# --- CPLEX cannot handle this, maybe need manual benders decomposition scheme.

function build_block_lopf(network, solver; formulation::String="angles_linear",
                         investment_type::String="continuous", decomposition::String="benders")

    m = build_lopf(network, solver;
                   formulation=formulation,
                   investment_type=investment_type, blockmodel=true)

    if decomposition == "benders"

        println("Apply Benders Decomposition")

        function b_decomp(varname::Symbol, varid::Tuple)
            
            master_vars = Set([:LN_s_nom, :LN_inv, :G_p_nom, :LN_opt, :LK_p_nom, :SU_p_nom, :ST_e_nom])

            if in(varname, master_vars)
                return (:B_MASTER, 0)
            else
                return (:B_SP, 0)
            end

        end

        add_Benders_decomposition(m, b_decomp)

    elseif decomposition == "dantzig-wolfe"

        # TODO: not implemented yet: evaluate suitability + assign constraints to subproblems

        println("Apply Dantzig Wolfe Decomposition")

        function dw_decomp(constr_name, constr_id)
            if constr_name == :mc
                return (:DW_MASTER, 0)
            else
                return (:DW_SP, constr_id[1])
            end
        end

        add_Dantzig_Wolfe_decomposition(m, dw_decomp)

    else

        error("not implemented")
    
    end    

    return m

end