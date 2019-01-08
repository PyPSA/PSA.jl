using BlockDecomposition

"""Builds a linear optimal power flow model to use with BlockDecomposition.jl"""
function build_block_lopf(network, solver; formulation::String="angles_linear",
                         investment_type::String="integer", decomposition::String="benders", n_sub::Int64=1)

    m = build_lopf(network, solver;
                   formulation=formulation,
                   investment_type=investment_type, blockmodel=true)

    if decomposition == "benders"

        println("Apply Benders Decomposition")

        function b_decomp(varname::Symbol, varid::Tuple)

            master_vars = Set([:LN_s_nom, :LN_opt, :LN_inv, :G_p_nom, :LK_p_nom, :SU_p_nom, :ST_e_nom])

            if in(varname, master_vars)
                return (:B_MASTER, 0)
            else
                return (:B_SP, (varid[2])%n_sub)
            end

        end

        add_Benders_decomposition(m, b_decomp)

    elseif decomposition == "dantzig-wolfe"

        # TODO: implement

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

        error("$decomposition is not implemented!")
    
    end    

    return m

end