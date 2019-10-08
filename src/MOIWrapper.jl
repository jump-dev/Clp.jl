
using LinQuadOptInterface
using .ClpCInterface

import LinearAlgebra

const LQOI = LinQuadOptInterface
const MOI  = LQOI.MOI

const SUPPORTED_OBJECTIVES = [
    LQOI.Linear,
    LQOI.SinVar
]

const SUPPORTED_CONSTRAINTS = [
    (LQOI.Linear, LQOI.EQ),
    (LQOI.Linear, LQOI.LE),
    (LQOI.Linear, LQOI.GE),
    (LQOI.SinVar, LQOI.EQ),
    (LQOI.SinVar, LQOI.LE),
    (LQOI.SinVar, LQOI.GE),
    (LQOI.SinVar, LQOI.IV),
    (LQOI.VecVar, MOI.Nonnegatives),
    (LQOI.VecVar, MOI.Nonpositives),
    (LQOI.VecVar, MOI.Zeros),
    (LQOI.VecLin, MOI.Nonnegatives),
    (LQOI.VecLin, MOI.Nonpositives),
    (LQOI.VecLin, MOI.Zeros)
]

mutable struct Optimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    params::Dict{Symbol,Any}
    Optimizer(::Nothing) = new()
end

### Options

# map option name to C function
const optionmap = Dict(
   :PrimalTolerance => set_primal_tolerance,
   :DualTolerance => set_dual_tolerance,
   :DualObjectiveLimit => set_dual_objective_limit,
   :MaximumIterations => set_maximum_iterations,
   :MaximumSeconds => set_maximum_seconds,
   :LogLevel => set_log_level,
   :Scaling => scaling,
   :Perturbation => set_perturbation,
   )
# These options are set by using the ClpSolve object
const solveoptionmap = Dict(
   :PresolveType => set_presolve_type,
   :SolveType => set_solve_type,
   :InfeasibleReturn => set_infeasible_return,
   )

function Optimizer(;kwargs...)
    optimizer = Optimizer(nothing)
    optimizer.params = Dict{String,Any}()
    MOI.empty!(optimizer)
    for (name,value) in kwargs
        optimizer.params[Symbol(name)] = value
    end
    return optimizer
end

LQOI.LinearQuadraticModel(::Type{Optimizer},env) = ClpModel()

LQOI.supported_constraints(optimizer::Optimizer) = SUPPORTED_CONSTRAINTS
LQOI.supported_objectives(optimizer::Optimizer) = SUPPORTED_OBJECTIVES

"""
    replace_inf(x::Real)

Return `Inf` if `x>1e20`, `-Inf` if `x<-1e20`, and `x` otherwise.
"""
function replace_inf(x::Real)
    if x > 1e20
        return Inf
    elseif x < -1e20
        return -Inf
    else
        return x
    end
end


function LQOI.get_objectivesense(model::Optimizer)
    s = get_obj_sense(model.inner)
    if s == 1.0
        return MOI.MIN_SENSE
    elseif s == -1.0
        return MOI.MAX_SENSE
    else
        @assert iszero(s)
        return MOI.FEASIBILITY_SENSE
    end
end


function LQOI.change_variable_bounds!(instance::Optimizer, cols::Vector{Int},
        values::Vector{Float64}, senses::Vector)
    upperbounds = get_col_upper(instance.inner)
    lowerbounds = get_col_lower(instance.inner)
    for (col, value, sense) in zip(cols, values, senses)
        if sense == Cchar('U')
            upperbounds[col] = value
        elseif sense == Cchar('L')
            lowerbounds[col] = value
        else
            error("sense is Cchar('$(Char(sense))'), but only Cchar('U') " *
                  "Cchar('L') are supported.")
        end
    end
    chg_column_upper(instance.inner, upperbounds)
    chg_column_lower(instance.inner, lowerbounds)
end

function LQOI.get_variable_lowerbound(instance::Optimizer, col::Int)
    lower = get_col_lower(instance.inner)
    return replace_inf(lower[col])
end

function LQOI.get_variable_upperbound(instance::Optimizer, col::Int)
    upper = get_col_upper(instance.inner)
    return replace_inf(upper[col])
end

function LQOI.get_number_linear_constraints(instance::Optimizer)
    return get_num_rows(instance.inner)
end

"""
    append_row(instance::Optimizer, row::Int, lower::Float64, upper::Float64,
               rows::Vector{Int}, cols::Vector{Int}, coefs::Vector{Float64})

Given a sparse matrix in the triplet-form `(rows, cols, coefs)`, add row `row`
with upper bound `upper` and lower bound  `lower` to the instance `instance`.
"""
function append_row(instance::Optimizer, row::Int, lower::Float64,
                    upper::Float64, rows::Vector{Int}, cols::Vector{Int},
                    coefs::Vector{Float64})
    indices = if row == length(rows)
        rows[row]:length(cols)
    else
        rows[row]:(rows[row+1]-1)
    end
    add_row(instance.inner, Cint(length(indices)), Cint.(cols[indices] .- 1),
            coefs[indices], lower, upper)
end

function LQOI.add_linear_constraints!(instance::Optimizer, A::LQOI.CSRMatrix{Float64},
        senses::Vector{Cchar}, right_hand_sides::Vector{Float64})
    rows = A.row_pointers
    cols = A.columns
    coefs = A.coefficients
    for (row, (rhs, sense)) in enumerate(zip(right_hand_sides, senses))
        if rhs > 1e20
            error("rhs must always be less than 1e20")
        elseif rhs < -1e20
            error("rhs must always be greater than -1e20")
        end
        lower = -Inf
        upper = Inf
        if sense == Cchar('L')
            upper = rhs
        elseif sense == Cchar('G')
            lower = rhs
        elseif sense == Cchar('E')
            upper = lower = rhs
        else
            error("sense must be Cchar(x) where x is in ['L','G',E']")
        end
        append_row(instance, row, lower, upper, rows, cols, coefs)
    end
end

function LQOI.add_ranged_constraints!(instance::Clp.Optimizer,
        A::LinQuadOptInterface.CSRMatrix{Float64},
        lb::Vector{Float64}, ub::Vector{Float64})
    rows = A.row_pointers
    cols = A.columns
    coefs = A.coefficients
    for row in 1:length(lb)
        append_row(instance, row, lb[row], ub[row], rows, cols, coefs)
    end
end

function LQOI.get_rhs(instance::Optimizer, row::Int)
    lower_bounds = get_row_lower(instance.inner)
    upper_bounds = get_row_upper(instance.inner)
    lower_bound = replace_inf(lower_bounds[row])
    upper_bound = replace_inf(upper_bounds[row])
    if lower_bound > -Inf
        return lower_bound
    elseif upper_bound < Inf
        return upper_bound
    else
        error("Either row_lower or row_upper must be of abs less than 1e20")
    end
end

function LQOI.get_linear_constraint(instance::Optimizer, row::Int)::Tuple{Vector{Int}, Vector{Float64}}
    A = get_constraint_matrix(instance.inner)
    A_row = A[row,:]
    return Array{Int}(A_row.nzind), A_row.nzval
end

function LQOI.change_objective_coefficient!(instance::Optimizer, col, coef)
    objcoefs = get_obj_coefficients(instance.inner)
    objcoefs[col] = coef
    chg_obj_coefficients(instance.inner, objcoefs)
end

function LQOI.change_rhs_coefficient!(instance::Optimizer, row, coef)
    lower_bounds = get_row_lower(instance.inner)
    upper_bounds = get_row_upper(instance.inner)
    lower_bound = replace_inf(lower_bounds[row])
    upper_bound = replace_inf(upper_bounds[row])
    if lower_bound > -Inf
        lower_bounds[row] = coef
        chg_row_lower(instance.inner, lower_bounds)
    elseif upper_bound < Inf
        upper_bounds[row] = coef
        chg_row_upper(instance.inner, upper_bounds)
    else
        error("Either row_lower or row_upper must be of abs less than 1e20")
    end
end

function LQOI.delete_linear_constraints!(instance::Optimizer, start_row::Int, end_row::Int)
    delete_rows(instance.inner, [Cint(i-1) for i in start_row:end_row])
end

function LQOI.change_linear_constraint_sense!(instance::Optimizer, rows::Vector{Int}, senses::Vector{Cchar})
    lower = replace_inf.(get_row_lower(instance.inner))
    upper = replace_inf.(get_row_upper(instance.inner))
    for (sense, row) in zip(senses, rows)
        lb = lower[row]
        ub = upper[row]
        if lb > -Inf
            rhs = lb
        elseif ub < Inf
            rhs = ub
        else
            error("Either row_lower or row_upper must be of abs less than 1e20")
        end
        if sense == Cchar('G')
            lower[row] = rhs
            upper[row] = Inf
        elseif sense == Cchar('L')
            lower[row] = -Inf
            upper[row] = rhs
        elseif sense == Cchar('E')
            lower[row] = rhs
            upper[row] = rhs
        end
    end
    chg_row_upper(instance.inner, upper)
    chg_row_lower(instance.inner, lower)
end

function LQOI.set_linear_objective!(instance::Optimizer, cols::Vector{Int}, coefs::Vector{Float64})
    objective_coefficients = zeros(Float64, get_num_cols(instance.inner))
    for (col, coef) in zip(cols, coefs)
        objective_coefficients[col] += coef
    end
    chg_obj_coefficients(instance.inner, objective_coefficients)
end

function LQOI.change_objective_sense!(instance::Optimizer, sense::Symbol)
    if sense == :min
        set_obj_sense(instance.inner, 1.0)
    elseif sense == :max
        set_obj_sense(instance.inner, -1.0)
    else
        error("sense must be either :min or :max")
    end
end

function LQOI.get_linear_objective!(instance::Optimizer, x::Vector{Float64})
    copyto!(x, get_obj_coefficients(instance.inner))
end

function LQOI.solve_linear_problem!(instance::Optimizer)
    solveroptions = ClpSolve()
    model = instance.inner
    for (name, value) in instance.params
        if haskey(optionmap, name)
            optionmap[name](model, value)
        elseif haskey(solveoptionmap, name)
            solveoptionmap[name](solveroptions,value)
        else
            error("Unrecognized option: $name")
        end
    end
    initial_solve_with_options(instance.inner, solveroptions)
end

function LQOI.get_variable_primal_solution!(instance::Optimizer, x::Vector{Float64})
    copyto!(x, primal_column_solution(instance.inner))
end

function LQOI.get_linear_primal_solution!(instance::Optimizer, x::Vector{Float64})
    copyto!(x, primal_row_solution(instance.inner))
end

function LQOI.get_variable_dual_solution!(instance::Optimizer, x::Vector{Float64})
    copyto!(x, dual_column_solution(instance.inner))
end

function LQOI.get_linear_dual_solution!(instance::Optimizer, x::Vector{Float64})
    copyto!(x, dual_row_solution(instance.inner))
end

function LQOI.get_objective_value(instance::Optimizer)
    return objective_value(instance.inner)
end

function LQOI.get_farkas_dual!(instance::Optimizer, result::Vector{Float64})
    copyto!(result, infeasibility_ray(instance.inner))
    LinearAlgebra.rmul!(result, -1.0)
end

function LQOI.get_unbounded_ray!(instance::Optimizer, result::Vector{Float64})
    copyto!(result, unbounded_ray(instance.inner))
end

function LQOI.get_termination_status(instance::Optimizer)
    status = ClpCInterface.status(instance.inner)
    if status == 0
        return MOI.OPTIMAL
    elseif status == 1
        return MOI.INFEASIBLE
    elseif status == 2
        return MOI.DUAL_INFEASIBLE
    elseif status == 3
        return MOI.OTHER_LIMIT
    elseif status == 4
        return MOI.OTHER_ERROR
    else
        error("status returned was $(status), but it must be in [0,1,2,3,4]")
    end
end

function LQOI.get_primal_status(instance::Optimizer)
    if is_proven_dual_infeasible(instance.inner)
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif primal_feasible(instance.inner)
        return MOI.FEASIBLE_POINT
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

function LQOI.get_dual_status(instance::Optimizer)
    if is_proven_primal_infeasible(instance.inner)
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif dual_feasible(instance.inner)
        return MOI.FEASIBLE_POINT
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

function LQOI.get_number_variables(instance::Optimizer)
    return get_num_cols(instance.inner)
end

function LQOI.add_variables!(instance::Optimizer, number_of_variables::Int)
    for i in 1:number_of_variables
        add_column(instance.inner, Cint(0), Int32[], Float64[], -Inf, Inf, 0.0)
    end
end

function LQOI.delete_variables!(instance::Optimizer, start_col::Int, end_col::Int)
    delete_columns(instance.inner, [Cint(i-1) for i in start_col:end_col])
end

# Corresponds to the `Status` struct defined in https://github.com/coin-or/Clp/blob/8419e63/Clp/src/ClpSimplex.hpp#L114.
const STATMAP = Dict(0x00 => MOI.BASIC, 0x01 => MOI.BASIC, 0x02 => MOI.NONBASIC_AT_UPPER,
                    0x03 => MOI.NONBASIC_AT_LOWER, 0x04 => MOI.SUPER_BASIC, 0x05 => MOI.NONBASIC)

function MOI.get(instance::Optimizer, ::MOI.ConstraintBasisStatus,
        ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}) where S <: Union{LQOI.LE, LQOI.GE, LQOI.EQ}
    row = instance[ci]
    stat = STATMAP[get_row_status(instance.inner, row)]
    # Single sided constraints should not specify `NONBASIC_AT_X` but only `NONBASIC`.
    if stat == MOI.NONBASIC_AT_LOWER || stat == MOI.NONBASIC_AT_UPPER
        return MOI.NONBASIC
    end
    return stat
end

function MOI.get(instance::Optimizer, ::MOI.ConstraintBasisStatus,
        vi::MOI.ConstraintIndex{MOI.SingleVariable, S}) where S <: Union{LQOI.LE, LQOI.GE, LQOI.EQ, LQOI.IV}
    col = instance.variable_mapping[instance[vi]]
    stat = STATMAP[get_column_status(instance.inner, col)]
    # If, e.g., a column is `NONBASIC_AT_LOWER` then the ≤ constraint is `BASIC`.
    if (S <: LQOI.LE && stat == MOI.NONBASIC_AT_LOWER) || (S <: LQOI.GE && stat == MOI.NONBASIC_AT_UPPER)
        return MOI.BASIC
    end
    # Single sided constraints should not specify `NONBASIC_AT_X` but only `NONBASIC`.
    if !(S <: LQOI.IV) && (stat == MOI.NONBASIC_AT_LOWER || stat == MOI.NONBASIC_AT_UPPER)
        return MOI.NONBASIC
    end
    return stat
end
