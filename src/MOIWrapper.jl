import MathOptInterface
import SparseArrays
using .ClpCInterface

const MOI = MathOptInterface

# Supported scalar sets
const SCALAR_SETS = Union{
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64}
}

# Maps Clp's parameters to getter/setter function
const CLP_OPTION_MAP = Dict(
    :PrimalTolerance => (Clp.primal_tolerance, Clp.set_primal_tolerance),
    :DualTolerance => (Clp.dual_tolerance, Clp.set_dual_tolerance),
    :DualObjectiveLimit => (Clp.dual_objective_limit, Clp.set_dual_objective_limit),
    :MaximumIterations => (Clp.maximum_iterations, Clp.set_maximum_iterations),
    :MaximumSeconds => (Clp.maximum_seconds, Clp.set_maximum_seconds),
    :LogLevel => (Clp.log_level, Clp.set_log_level),
    :Scaling => (Clp.scaling_flag, Clp.scaling),
    :Perturbation => (Clp.perturbation, Clp.set_perturbation),
    :Algorithm => (Clp.algorithm, Clp.set_algorithm)
)

const SOLVE_OPTION_MAP = Dict(
   :PresolveType => (Clp.get_presolve_type, Clp.set_presolve_type),
   :SolveType => (Clp.get_solve_type, Clp.set_solve_type),
   :InfeasibleReturn => (Clp.infeasible_return, Clp.set_infeasible_return)
)

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Clp.ClpModel
    solver_options::Clp.ClpSolve
    options_set::Set{Symbol}
    optimize_called::Bool
    function Optimizer(;kwargs...)
        inner_model = Clp.ClpModel()
        solver_options = Clp.ClpSolve()
        model = new(inner_model, solver_options, Set{Symbol}(), false)
        for (key, value) in kwargs
            MOI.set(model, MOI.RawParameter(String(key)), value)
        end
        return model
    end
end

# ====================
#   empty functions
# ====================

function MOI.is_empty(model::Optimizer)
    # A problem is empty if it has no variable and no linear constraints
    return (Clp.get_num_rows(model.inner) == 0) && (Clp.get_num_cols(model.inner) == 0)
end

function MOI.empty!(model::Optimizer)
    old_model = model.inner
    model.inner = Clp.ClpModel()
    # Copy parameters from old model into new model
    for key in model.options_set
        getter, setter = CLP_OPTION_MAP[key]
        setter(model.inner, getter(old_model))
    end
    model.optimize_called = false
    return
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Clp"

# TODO: improve type-stability of `MOI.RawParameter`-related methods.
MOI.supports(::Optimizer, param::MOI.RawParameter) = true

function MOI.set(model::Optimizer, param::MOI.RawParameter, value)
    key = Symbol(param.name)
    if haskey(CLP_OPTION_MAP, key)
        push!(model.options_set, key)
        CLP_OPTION_MAP[key][2](model.inner, value)
    elseif haskey(SOLVE_OPTION_MAP, key)
        SOLVE_OPTION_MAP[key][2](model.solver_options, value)
    else
        throw(MOI.UnsupportedAttribute(param))
    end
    return
end

function MOI.get(model::Optimizer, param::MOI.RawParameter)
    key = Symbol(param.name)
    if haskey(CLP_OPTION_MAP, key)
        return CLP_OPTION_MAP[key][1](model.inner)
    elseif haskey(SOLVE_OPTION_MAP, key)
        return SOLVE_OPTION_MAP[key][1](model.solver_options)
    else
        throw(MOI.UnsupportedAttribute(param))
    end
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    push!(model.options_set, :LogLevel)
    if value
        Clp.set_log_level(model.inner, 0)
    else
        Clp.set_log_level(model.inner, 1)
    end
    return
end

function MOI.get(model::Optimizer, ::MOI.Silent)
    return Clp.log_level(model.inner) == 0
end

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value)
    push!(model.options_set, :MaximumSeconds)
    if value === nothing
        # De-activate time limit
        Clp.set_maximum_seconds(model.inner, -1.0)
    else
        Clp.set_maximum_seconds(model.inner, value)
    end
    return
end

# Clp behaves weird here
MOI.get(model::Optimizer, ::MOI.TimeLimitSec) = Clp.maximum_seconds(model.inner)

MOI.supports(::Optimizer, ::MOI.NumberOfThreads) = false


# ========================================
#   Supported constraints and objectives
# ========================================

function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{<:SCALAR_SETS}
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.SingleVariable}, ::Type{<:SCALAR_SETS}
)
    return true
end

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.supports(
    ::Optimizer, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}
)
    return true
end

# =======================
#   `copy_to` function
# =======================

_add_bounds(::Vector{Float64}, ub, i, s::MOI.LessThan{Float64}) = ub[i] = s.upper
_add_bounds(lb, ::Vector{Float64}, i, s::MOI.GreaterThan{Float64}) = lb[i] = s.lower
_add_bounds(lb, ub, i, s::MOI.EqualTo{Float64}) = lb[i], ub[i] = s.value, s.value
_add_bounds(lb, ub, i, s::MOI.Interval{Float64}) = lb[i], ub[i] = s.lower, s.upper

function _extract_bound_data(src, mapping, lb, ub, S)
    for con_index in MOI.get(
        src, MOI.ListOfConstraintIndices{MOI.SingleVariable, S}()
    )
        f = MOI.get(src, MOI.ConstraintFunction(), con_index)
        s = MOI.get(src, MOI.ConstraintSet(), con_index)
        column = mapping.varmap[f.variable].value
        _add_bounds(lb, ub, column, s)
        mapping.conmap[con_index] = MOI.ConstraintIndex{MOI.SingleVariable, S}(column)
    end
end

function _copy_to_columns(dest, src, mapping)
    x_src = MOI.get(src, MOI.ListOfVariableIndices())
    N = Cint(length(x_src))
    for i = 1:N
        mapping.varmap[x_src[i]] = MOI.VariableIndex(i)
    end

    fobj = MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    c = fill(0.0, N)
    for term in fobj.terms
        i = mapping.varmap[term.variable_index].value
        c[i] += term.coefficient
    end
    # Clp seems to negates the objective offset
    Clp.set_objective_offset(dest.inner, -fobj.constant)
    return N, c
end

_bounds(s::MOI.GreaterThan{Float64}) = (s.lower, Inf)
_bounds(s::MOI.LessThan{Float64}) = (-Inf, s.upper)
_bounds(s::MOI.EqualTo{Float64}) = (s.value, s.value)
_bounds(s::MOI.Interval{Float64}) = (s.lower, s.upper)

function _extract_row_data(src, mapping, lb, ub, I, J, V, S)
    row = length(I) == 0 ? 1 : I[end] + 1
    for c_index in MOI.get(
        src, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64}, S}()
    )
        f = MOI.get(src, MOI.ConstraintFunction(), c_index)
        l, u = _bounds(MOI.get(src, MOI.ConstraintSet(), c_index))
        push!(lb, l - f.constant)
        push!(ub, u - f.constant)
        for term in f.terms
            push!(I, row)
            push!(J, Cint(mapping.varmap[term.variable_index].value))
            push!(V, term.coefficient)
        end
        mapping.conmap[c_index] = MOI.ConstraintIndex{
            MOI.ScalarAffineFunction{Float64}, S
        }(length(ub))
        row += 1
    end
    return
end

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false
)
    @assert MOI.is_empty(dest)
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if !MOI.supports_constraint(dest, F, S)
            throw(MOI.UnsupportedConstraint{F, S}("Clp.Optimizer does not support constraints of type $F-in-$S."))
        end
    end
    fobj_type = MOI.get(src, MOI.ObjectiveFunctionType())
    if !MOI.supports(dest, MOI.ObjectiveFunction{fobj_type}())
        throw(MOI.UnsupportedAttribute(MOI.ObjectiveFunction(fobj_type)))
    end
    mapping = MOI.Utilities.IndexMap()
    N, c = _copy_to_columns(dest, src, mapping)
    cl, cu = fill(-Inf, N), fill(Inf, N)
    rl, ru, I, J, V = Float64[], Float64[], Cint[], Cint[], Float64[]
    for S in (
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.Interval{Float64},
    )
        _extract_bound_data(src, mapping, cl, cu, S)
        _extract_row_data(src, mapping, rl, ru, I, J, V, S)
    end
    M = Cint(length(rl))
    A = SparseArrays.sparse(I, J, V, M, N)
    Clp.load_problem(
        dest.inner,
        A.n,
        A.m,
        A.colptr .- Cint(1),
        A.rowval .- Cint(1),
        A.nzval,
        cl,
        cu,
        c,
        rl,
        ru
    )

    sense = MOI.get(src, MOI.ObjectiveSense())
    if sense == MOI.MIN_SENSE
        Clp.set_obj_sense(dest.inner, 1)
    elseif sense == MOI.MAX_SENSE
        Clp.set_obj_sense(dest.inner, -1)
    else
        @assert sense == MOI.FEASIBILITY_SENSE
        Clp.set_obj_sense(dest.inner, 0)
    end
    return mapping
end

# ===============================
#   Optimize and post-optimize
# ===============================

function MOI.optimize!(model::Optimizer)
    Clp.initial_solve_with_options(model.inner, model.solver_options)
    model.optimize_called = true
    return
end

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)
    return Clp.get_num_cols(model.inner)
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return Clp.get_obj_value(model.inner)
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if !model.optimize_called
        return MOI.OPTIMIZE_NOT_CALLED
    end
    st = Clp.status(model.inner)
    if st == 0
        return MOI.OPTIMAL
    elseif st == 1
        return MOI.INFEASIBLE
    elseif st == 2
        return MOI.DUAL_INFEASIBLE
    elseif st == 3
        # No more granular information that "some limit is reached"
        return MOI.OTHER_LIMIT
    else
        return MOI.OTHER_ERROR
    end
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    if !model.optimize_called
        return "MOI.OPTIMIZE_NOT_CALLED"
    end
    st = Clp.status(model.inner)
    if st == 0
        return "0 - optimal"
    elseif st == 1
        return "1 - primal infeasible"
    elseif st == 2
        return "2 - dual infeasible"
    elseif st == 3
        return "3 - stopped on iterations etc"
    elseif st == 4
        return "4 - stopped due to errors"
    else
        error("Expected integer in [0, 4] but got $st")
    end
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    if Clp.primal_feasible(model.inner)
        return 1
    elseif Clp.dual_feasible(model.inner)
        return 1
    elseif Clp.is_proven_primal_infeasible(model.inner)
        return 1
    elseif Clp.is_proven_dual_infeasible(model.inner)
        return 1
    end
    return 0
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif Clp.is_proven_dual_infeasible(model.inner)
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif Clp.primal_feasible(model.inner)
        return MOI.FEASIBLE_POINT
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif Clp.is_proven_primal_infeasible(model.inner)
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif Clp.dual_feasible(model.inner)
        return MOI.FEASIBLE_POINT
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

# ===================
#   Primal solution
# ===================

function MOI.get(
    model::Optimizer, attr::MOI.VariablePrimal, x::MOI.VariableIndex
)
    MOI.check_result_index_bounds(model, attr)
    primal_status = MOI.get(model, MOI.PrimalStatus())
    if primal_status == MOI.INFEASIBILITY_CERTIFICATE
        sol = Clp.unbounded_ray(model.inner)
        return sol[x.value]
    elseif primal_status == MOI.FEASIBLE_POINT
        sol = Clp.get_col_solution(model.inner)
        return sol[x.value]
    else
        error("Primal solution not available")
    end
end

function MOI.get(
    model::Optimizer, attr::MOI.VariablePrimal, xs::Vector{MOI.VariableIndex}
)
    MOI.check_result_index_bounds(model, attr)
    col_indices = [idx.value for idx in xs]
    primal_status = MOI.get(model, MOI.PrimalStatus())
    if primal_status == MOI.INFEASIBILITY_CERTIFICATE
        sol = Clp.unbounded_ray(model.inner)
        return sol[col_indices]
    elseif primal_status == MOI.FEASIBLE_POINT
        sol = Clp.get_col_solution(model.inner)
        return sol[col_indices]
    else
        error("Primal solution not available")
    end
end

# TODO: What happens if model is unbounded / infeasible?
function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:SCALAR_SETS}
)
    MOI.check_result_index_bounds(model, attr)
    sol = Clp.get_row_activity(model.inner)
    return sol[c.value]
end

# TODO: What happens if model is unbounded / infeasible?
function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.SingleVariable, <:SCALAR_SETS}
)
    MOI.check_result_index_bounds(model, attr)
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end


# =================
#   Dual solution
# =================

# If sense is maximize, we negate all the duals to follow MOI conventions
# Feasibility problems are treated as a minimization

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, <:SCALAR_SETS}
)
    MOI.check_result_index_bounds(model, attr)
    sense = (Clp.get_obj_sense(model.inner) == -1) ? -1 : 1
    dual_status = MOI.get(model, MOI.DualStatus())
    if dual_status == MOI.FEASIBLE_POINT
        dsol = Clp.get_row_price(model.inner)
        return sense * dsol[c.value]
    elseif dual_status == MOI.INFEASIBILITY_CERTIFICATE
        dsol = Clp.infeasibility_ray(model.inner)
        return sense * dsol[c.value]
    else
        error("Dual solution not available")
    end
end

# TODO: what happens if problem is unbounded / infeasible?
function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}
)
    MOI.check_result_index_bounds(model, attr)
    rc = Clp.get_reduced_cost(model.inner)[c.value]
    sense = (Clp.get_obj_sense(model.inner) == -1) ? -1 : 1
    if sense == 1 && rc <= 0.0
        return rc
    elseif sense == -1 && rc >= 0.0
        return -rc
    else
        return 0.0
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}
)
    MOI.check_result_index_bounds(model, attr)
    rc = Clp.get_reduced_cost(model.inner)[c.value]
    sense = (Clp.get_obj_sense(model.inner) == -1) ? -1 : 1
    if sense == 1 && rc >= 0.0
        return rc
    elseif sense == -1 && rc <= 0.0
        return -rc
    else
        return 0.0
    end
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{
        MOI.SingleVariable,
        <:Union{MOI.Interval{Float64}, MOI.EqualTo{Float64}}
    }
)
    MOI.check_result_index_bounds(model, attr)
    sense = (Clp.get_obj_sense(model.inner) == -1) ? -1 : 1
    return sense * Clp.get_reduced_cost(model.inner)[c.value]
end
