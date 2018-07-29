export ClpOptimizer

using LinQuadOptInterface
using .ClpCInterface

const LQOI = LinQuadOptInterface
const MOI  = LQOI.MOI

const SUPPORTED_OBJECTIVES = [
    LQOI.Linear
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

mutable struct ClpOptimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    params::Dict{Symbol,Any}
    ClpOptimizer(::Void) = new()
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

function setoption(optimizer::ClpOptimizer, name::Symbol, value)
    if haskey(optionmap, name)
        optionmap[name](optimizer.inner,value)
    elseif haskey(solveoptionmap, name)
        solveoptionmap[name](optimizer.solveroptions,value)
    else
        error("Unrecognized option: $name")
    end
end

function ClpOptimizer(;kwargs...)
    optimizer = ClpOptimizer(nothing)
    optimizer.params = Dict{String,Any}()
    MOI.empty!(optimizer)
    for (name,value) in kwargs
        optimizer.params[Symbol(name)] = value
    end
    return optimizer
end

LQOI.LinearQuadraticModel(::Type{ClpOptimizer},env) = ClpModel()

LQOI.supported_constraints(optimizer::ClpOptimizer) = SUPPORTED_CONSTRAINTS
LQOI.supported_objectives(optimizer::ClpOptimizer) = SUPPORTED_OBJECTIVES

"""
    replace_inf(x::Vector{<:Real})

Replaces any occurances of an element `x[i]>1e20` with `Inf` and `x[i]<-1e20`
with `-Inf`.
"""
function replace_inf(x::Vector{<:Real})
    for i in 1:length(x)
        if x[i] > 1e20
            x[i] = Inf
        elseif x[i] < -1e20
            x[i] = -Inf
        end
    end
    return x
end

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

function LQOI.change_variable_bounds!(instance::ClpOptimizer, cols::Vector{Int},
        values::Vector{Float64}, senses::Vector)
    upperbounds = replace_inf(get_col_upper(instance.inner))
    lowerbounds = replace_inf(get_col_lower(instance.inner))
    for i in 1:length(senses)
        if senses[i] == Cchar('U')
            upperbounds[cols[i]] = values[i]
        elseif senses[i] == Cchar('L')
            lowerbounds[cols[i]] = values[i]
        else
            error("sense is Cchar('$(Char(senses[i]))'), but only Cchar('U') " *
                  "Cchar('L') are supported.")
        end
    end
    chg_column_upper(instance.inner, upperbounds)
    chg_column_lower(instance.inner, lowerbounds)
end

function LQOI.get_variable_lowerbound(instance::ClpOptimizer, col::Int)
    lower = get_col_lower(instance.inner)
    return replace_inf(lower[col])
end

function LQOI.get_variable_upperbound(instance::ClpOptimizer, col::Int)
    upper = get_col_upper(instance.inner)
    return replace_inf(upper[col])
end

function LQOI.get_number_linear_constraints(instance::ClpOptimizer)
    return get_num_rows(instance.inner)
end

"""
    append_row(instance::ClpOptimizer, row::Int, lower::Float64, upper::Float64,
               rows::Vector{Int}, cols::Vector{Int}, coefs::Vector{Float64})

Given a sparse matrix in the triplet-form `(rows, cols, coefs)`, add row `row`
with upper bound `upper` and lower bound  `lower` to the instance `instance`.
"""
function append_row(instance::ClpOptimizer, row::Int, lower::Float64,
                    upper::Float64, rows::Vector{Int}, cols::Vector{Int},
                    coefs::Vector{Float64})
    row_variables = Int32[]
    row_coefficients = Float64[]
    first = rows[row]
    last = (row==length(rows)) ? length(cols) : rows[row+1]-1
    for i in first:last
        push!(row_coefficients, coefs[i])
        push!(row_variables, cols[i] - 1)
    end
    add_row(instance.inner, Cint(length(row_variables)), row_variables,
            row_coefficients, lower, upper)
end

function LQOI.add_linear_constraints!(instance::ClpOptimizer, A::LQOI.CSRMatrix{Float64},
        senses::Vector{Cchar}, right_hand_sides::Vector{Float64})
    rows = A.row_pointers
    cols = A.columns
    coefs = A.coefficients
    for (row, (rhs, sense)) in enumerate(zip(right_hand_sides, senses))
        if (rhs > 1e20)
            error("rhs must always be less than 1e20")
        elseif (rhs < -1e20)
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

function LQOI.add_ranged_constraints!(instance::Clp.ClpOptimizer,
        A::LinQuadOptInterface.CSRMatrix{Float64},
        lb::Vector{Float64}, ub::Vector{Float64})
    rows = A.row_pointers
    cols = A.columns
    coefs = A.coefficients
    for row in 1:length(lb)
        append_row(instance, row, lb[row], ub[row], rows, cols, coefs)
    end
end

function LQOI.get_rhs(instance::ClpOptimizer, row::Int)
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

function LQOI.get_linear_constraint(instance::ClpOptimizer, row::Int)::Tuple{Vector{Int}, Vector{Float64}}
    A = get_constraint_matrix(instance.inner)
    A_row = A[row,:]
    return (Array{Int}(A_row.nzind), A_row.nzval)
end

function LQOI.change_objective_coefficient!(instance::ClpOptimizer, col, coef)
    objcoefs = get_obj_coefficients(instance.inner)
    objcoefs[col] = coef
    chg_obj_coefficients(instance.inner, objcoefs)
end

function LQOI.change_rhs_coefficient!(instance::ClpOptimizer, row, coef)
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

function LQOI.delete_linear_constraints!(instance::ClpOptimizer, start_row::Int, end_row::Int)
    rows = [Cint(i-1) for i in start_row:end_row]
    delete_rows(instance.inner, rows)
end

function LQOI.change_linear_constraint_sense!(instance::ClpOptimizer, rows::Vector{Int}, senses::Vector{Cchar})
    lower = replace_inf(get_row_lower(instance.inner))
    upper = replace_inf(get_row_upper(instance.inner))
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
            lb = rhs
            ub = Inf
        elseif sense == Cchar('L')
            lb = -Inf
            ub = rhs
        elseif sense == Cchar('E')
            lb = rhs
            ub = rhs
        end
        lower[row] = lb
        upper[row] = ub
    end
    chg_row_upper(instance.inner, upper)
    chg_row_lower(instance.inner, lower)
end

function LQOI.set_linear_objective!(instance::ClpOptimizer, cols::Vector{Int}, coefs::Vector{Float64})
    objective_coefficients = zeros(Float64, get_num_cols(instance.inner))
    for (col, coef) in zip(cols, coefs)
        objective_coefficients[col] += coef
    end
    chg_obj_coefficients(instance.inner, objective_coefficients)
end

function LQOI.change_objective_sense!(instance::ClpOptimizer, sense::Symbol)
    if sense == :min
        set_obj_sense(instance.inner, 1.0)
    elseif sense == :max
        set_obj_sense(instance.inner, -1.0)
    else
        error("sense must be either :min or :max")
    end
end

function LQOI.get_linear_objective!(instance::ClpOptimizer, x::Vector{Float64})
    copy!(x, get_obj_coefficients(instance.inner))
end

function LQOI.solve_linear_problem!(instance::ClpOptimizer)
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

function LQOI.get_variable_primal_solution!(instance::ClpOptimizer, x::Vector{Float64})
    copy!(x, primal_column_solution(instance.inner))
end

function LQOI.get_linear_primal_solution!(instance::ClpOptimizer, x::Vector{Float64})
    copy!(x, primal_row_solution(instance.inner))
end

function LQOI.get_variable_dual_solution!(instance::ClpOptimizer, x::Vector{Float64})
    copy!(x, dual_column_solution(instance.inner))
end

function LQOI.get_linear_dual_solution!(instance::ClpOptimizer, x::Vector{Float64})
    copy!(x, dual_row_solution(instance.inner))
end

function LQOI.get_objective_value(instance::ClpOptimizer)
    return objective_value(instance.inner)
end

function LQOI.get_farkas_dual!(instance::ClpOptimizer, result::Vector{Float64})
    copy!(result, infeasibility_ray(instance.inner))
    scale!(result, -1.0)
end

function LQOI.get_unbounded_ray!(instance::ClpOptimizer, result::Vector{Float64})
    copy!(result, unbounded_ray(instance.inner))
end

function LQOI.get_termination_status(instance::ClpOptimizer)
    status = ClpCInterface.status(instance.inner)
    if status == 0
        return MOI.Success
    elseif status == 1
        if is_proven_primal_infeasible(instance.inner)
            return MOI.Success
        else
            return MOI.InfeasibleNoResult
        end
    elseif status == 2
        if is_proven_dual_infeasible(instance.inner)
            return MOI.Success
        else
            return MOI.UnboundedNoResult
        end
    elseif status == 3
        return MOI.OtherLimit
    elseif status == 4
        return MOI.OtherError
    else
        error("status returned was $(status), but it must be in [0,1,2,3,4]")
    end
end

function LQOI.get_primal_status(instance::ClpOptimizer)
    if is_proven_dual_infeasible(instance.inner)
        return MOI.InfeasibilityCertificate
    elseif primal_feasible(instance.inner)
        return MOI.FeasiblePoint
    else
        return MOI.UnknownResultStatus
    end
end

function LQOI.get_dual_status(instance::ClpOptimizer)
    if is_proven_primal_infeasible(instance.inner)
        return MOI.InfeasibilityCertificate
    elseif dual_feasible(instance.inner)
        return MOI.FeasiblePoint
    else
        return MOI.UnknownResultStatus
    end
end

function LQOI.get_number_variables(instance::ClpOptimizer)
    return get_num_cols(instance.inner)
end

function LQOI.add_variables!(instance::ClpOptimizer, n::Int)
    for i in 1:n
        add_column(instance.inner, Cint(0), Int32[], Float64[], -Inf, Inf, 0.0)
    end
end

function LQOI.delete_variables!(instance::ClpOptimizer, start_col::Int, end_col::Int)
    columns = [Cint(i-1) for i in start_col:end_col]
    delete_columns(instance.inner, columns)
end
