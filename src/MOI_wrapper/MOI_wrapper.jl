import MathOptInterface
import SparseArrays

const MOI = MathOptInterface

# Supported scalar sets
const SCALAR_SETS = Union{
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64},
}

mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Ptr{Cvoid}
    solver_options::Ptr{Cvoid}
    options_set::Set{String}
    optimize_called::Bool
    solve_time::Float64
    # Work-around for upstream bug in Clp:
    maximumSeconds::Float64

    """
        Optimizer()

    Create a new Optimizer object.

    Set optimizer attributes using `MOI.RawParameter` or
    `JuMP.set_optimizer_atttribute`. For a list of supported parameter names,
    see `Clp.SUPPORTED_PARAMETERS`.

    ## Example

        using JuMP, Clp
        model = JuMP.Model(Clp.Optimizer)
        set_optimizer_attribute(model, "LogLevel", 0)
    """
    function Optimizer(; kwargs...)
        model = new(
            Clp_newModel(),
            ClpSolve_new(),
            Set{String}(),
            false,
            0.0,
            -1.0,
        )
        if length(kwargs) > 0
            @warn("""Passing optimizer attributes as keyword arguments to
            Clp.Optimizer is deprecated. Use
                MOI.set(model, MOI.RawParameter("key"), value)
            or
                JuMP.set_optimizer_attribute(model, "key", value)
            instead.
            """)
        end
        for (key, value) in kwargs
            MOI.set(model, MOI.RawParameter(String(key)), value)
        end
        finalizer(model) do m
            Clp_deleteModel(m.inner)
            ClpSolve_delete(m.solver_options)
        end
        return model
    end
end

# ====================
#   empty functions
# ====================

function MOI.is_empty(model::Optimizer)
    # A problem is empty if it has no variable and no linear constraints
    return (Clp_getNumRows(model.inner) == 0) && (Clp_getNumCols(model.inner) == 0)
end

function MOI.empty!(model::Optimizer)
    old_model = model.inner
    model.inner = Clp_newModel()
    model.optimize_called = false
    model.solve_time = 0.0
    # Copy parameters from old model into new model
    old_options = copy(model.options_set)
    empty!(model.options_set)
    for key in old_options
        value = MOI.get(model, MOI.RawParameter(key))
        MOI.set(model, MOI.RawParameter(key), value)
    end
    Clp_deleteModel(old_model)
    # Work-around for maximumSeconds
    Clp_setMaximumSeconds(model.inner, model.maximumSeconds)
    return
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Clp"

# If you update this list, remember to update the README documentation!
const SUPPORTED_PARAMETERS = (
    "PrimalTolerance",
    "DualTolerance",
    "DualObjectiveLimit",
    "MaximumIterations",
    "MaximumSeconds",
    "LogLevel",
    "Scaling",
    "Perturbation",
    # TODO(odow): `Algorithm` is excluded from the README because it isn't
    # apparent what it controls. Use `SolveType` instead.
    "Algorithm",
    "PresolveType",
    "SolveType",
    "InfeasibleReturn",
)

function MOI.supports(::Optimizer, param::MOI.RawParameter)
    return param.name in SUPPORTED_PARAMETERS
end

function MOI.set(model::Optimizer, param::MOI.RawParameter, value)
    name = String(param.name)
    push!(model.options_set, name)
    if name == "PrimalTolerance"
        Clp_setPrimalTolerance(model.inner, value)
    elseif name == "DualTolerance"
        Clp_setDualTolerance(model.inner, value)
    elseif name == "DualObjectiveLimit"
        Clp_setDualObjectiveLimit(model.inner, value)
    elseif name == "MaximumIterations"
        Clp_setMaximumIterations(model.inner, value)
    elseif name == "MaximumSeconds"
        Clp_setMaximumSeconds(model.inner, value)
    elseif name == "LogLevel"
        Clp_setLogLevel(model.inner, value)
    elseif name == "Scaling"
        Clp_scaling(model.inner, value)
    elseif name == "Perturbation"
        Clp_setPerturbation(model.inner, value)
    elseif name == "Algorithm"
        Clp_setAlgorithm(model.inner, value)
    elseif name == "PresolveType"
        ClpSolve_setPresolveType(model.solver_options, value, -1)
    elseif name == "SolveType"
        ClpSolve_setSolveType(model.solver_options, value, -1)
    elseif name == "InfeasibleReturn"
        ClpSolve_setInfeasibleReturn(model.solver_options, value)
    else
        throw(MOI.UnsupportedAttribute(param))
    end
    return
end

function MOI.get(model::Optimizer, param::MOI.RawParameter)
    name = String(param.name)
    if name == "PrimalTolerance"
        return Clp_primalTolerance(model.inner)
    elseif name == "DualTolerance"
        return Clp_dualTolerance(model.inner)
    elseif name == "DualObjectiveLimit"
        return Clp_dualObjectiveLimit(model.inner)
    elseif name == "MaximumIterations"
        # TODO(odow): Clp doesn't follow convention with maximumIterations!
        return maximumIterations(model.inner)
    elseif name == "MaximumSeconds"
        return Clp_maximumSeconds(model.inner)
    elseif name == "LogLevel"
        return Clp_logLevel(model.inner)
    elseif name == "Scaling"
        return Clp_scalingFlag(model.inner)
    elseif name == "Perturbation"
        return Clp_perturbation(model.inner)
    elseif name == "Algorithm"
        return Clp_algorithm(model.inner)
    elseif name == "PresolveType"
        return ClpSolve_getPresolveType(model.solver_options)
    elseif name == "SolveType"
        return ClpSolve_getSolveType(model.solver_options)
    elseif name == "InfeasibleReturn"
        return ClpSolve_infeasibleReturn(model.solver_options)
    else
        throw(MOI.UnsupportedAttribute(param))
    end
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    push!(model.options_set, "LogLevel")
    Clp_setLogLevel(model.inner, value ? 0 : 1)
    return
end

function MOI.get(model::Optimizer, ::MOI.Silent)
    return Clp_logLevel(model.inner) == 0
end

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value)
    push!(model.options_set, "MaximumSeconds")
    value = value === nothing ? -1.0 : value
    Clp_setMaximumSeconds(model.inner, value)
    model.maximumSeconds = value
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    # TODO(odow): replace with `Clp_maximumSeconds(model.inner)` when upstream
    # is fixed.
    return model.maximumSeconds
end

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

function _extract_bound_data(src, mapping, lb, ub, ::Type{S}) where S
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
    Clp_setObjectiveOffset(dest.inner, -fobj.constant)
    return N, c
end

_bounds(s::MOI.GreaterThan{Float64}) = (s.lower, Inf)
_bounds(s::MOI.LessThan{Float64}) = (-Inf, s.upper)
_bounds(s::MOI.EqualTo{Float64}) = (s.value, s.value)
_bounds(s::MOI.Interval{Float64}) = (s.lower, s.upper)

function add_sizehint!(vec, n)
    len = length(vec)
    return sizehint!(vec, len + n)
end

function _extract_row_data(src, mapping, lb, ub, I, J, V, ::Type{S}) where S
    row = length(I) == 0 ? 1 : I[end] + 1
    list = MOI.get(
        src, MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Float64}, S}()
    )
    add_sizehint!(lb, length(list))
    add_sizehint!(ub, length(list))
    n_terms = 0
    fs = Array{MOI.ScalarAffineFunction{Float64}}(undef, length(list))
    for (i,c_index) in enumerate(list)
        f = MOI.get(src, MOI.ConstraintFunction(), c_index)
        fs[i] = f
        l, u = _bounds(MOI.get(src, MOI.ConstraintSet(), c_index))
        push!(lb, l - f.constant)
        push!(ub, u - f.constant)
        n_terms += length(f.terms)
    end
    add_sizehint!(I, n_terms)
    add_sizehint!(J, n_terms)
    add_sizehint!(V, n_terms)
    for (i,c_index) in enumerate(list)
        f = fs[i]#MOI.get(src, MOI.ConstraintFunction(), c_index)
        for term in f.terms
            push!(I, row)
            push!(J, Cint(mapping.varmap[term.variable_index].value))
            push!(V, term.coefficient)
        end
        mapping.conmap[c_index] = MOI.ConstraintIndex{
            MOI.ScalarAffineFunction{Float64}, S
        }(row)
        row += 1
    end
    return
end

function test_data(src, dest)
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        if !MOI.supports_constraint(dest, F, S)
            throw(MOI.UnsupportedConstraint{F, S}("Clp.Optimizer does not support constraints of type $F-in-$S."))
        end
    end
    fobj_type = MOI.get(src, MOI.ObjectiveFunctionType())
    if !MOI.supports(dest, MOI.ObjectiveFunction{fobj_type}())
        throw(MOI.UnsupportedAttribute(MOI.ObjectiveFunction(fobj_type)))
    end
end

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.ModelLike;
    copy_names::Bool = false
)
    @assert MOI.is_empty(dest)
    test_data(src, dest)

    mapping = MOI.Utilities.IndexMap()
    N, c = _copy_to_columns(dest, src, mapping)
    cl, cu = fill(-Inf, N), fill(Inf, N)
    rl, ru, I, J, V = Float64[], Float64[], Cint[], Cint[], Float64[]

    _extract_bound_data(src, mapping, cl, cu, MOI.GreaterThan{Float64})
    _extract_row_data(src, mapping, rl, ru, I, J, V, MOI.GreaterThan{Float64})
    _extract_bound_data(src, mapping, cl, cu, MOI.LessThan{Float64})
    _extract_row_data(src, mapping, rl, ru, I, J, V, MOI.LessThan{Float64})
    _extract_bound_data(src, mapping, cl, cu, MOI.EqualTo{Float64})
    _extract_row_data(src, mapping, rl, ru, I, J, V, MOI.EqualTo{Float64})
    _extract_bound_data(src, mapping, cl, cu, MOI.Interval{Float64})
    _extract_row_data(src, mapping, rl, ru, I, J, V, MOI.Interval{Float64})

    M = Cint(length(rl))
    A = SparseArrays.sparse(I, J, V, M, N)
    Clp_loadProblem(
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
        Clp_setObjSense(dest.inner, 1)
    elseif sense == MOI.MAX_SENSE
        Clp_setObjSense(dest.inner, -1)
    else
        @assert sense == MOI.FEASIBILITY_SENSE
        Clp_setObjSense(dest.inner, 0)
    end
    return mapping
end

# ===============================
#   Optimize and post-optimize
# ===============================

function MOI.optimize!(model::Optimizer)
    t = time()
    Clp_initialSolveWithOptions(model.inner, model.solver_options)
    model.solve_time = time() - t
    model.optimize_called = true
    return
end

function MOI.get(model::Optimizer, ::MOI.SolveTime)
    return model.solve_time
end

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)
    return Clp_getNumCols(model.inner)
end

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    return Clp_getObjValue(model.inner)
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if !model.optimize_called
        return MOI.OPTIMIZE_NOT_CALLED
    end
    st = Clp_status(model.inner)
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
    st = Clp_status(model.inner)
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
    if Clp_primalFeasible(model.inner) != 0
        return 1
    elseif Clp_dualFeasible(model.inner) != 0
        return 1
    elseif Clp_isProvenPrimalInfeasible(model.inner) != 0
        return 1
    elseif Clp_isProvenDualInfeasible(model.inner) != 0
        return 1
    end
    return 0
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif Clp_isProvenDualInfeasible(model.inner) != 0
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif Clp_primalFeasible(model.inner) != 0
        return MOI.FEASIBLE_POINT
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if attr.N != 1
        return MOI.NO_SOLUTION
    elseif Clp_isProvenPrimalInfeasible(model.inner) != 0
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif Clp_dualFeasible(model.inner) != 0
        return MOI.FEASIBLE_POINT
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

# ===================
#   Primal solution
# ===================

function _unsafe_wrap_clp_array(
    model::Ptr{Cvoid},
    f::Function,
    n::Integer,
    indices;
    own::Bool = false
)
    p = f(model)
    if p == C_NULL
        return map(x -> NaN, indices)
    end
    x = unsafe_wrap(Array, p, (n,); own = own)
    return x[indices]
end

function MOI.get(
    model::Optimizer, attr::MOI.VariablePrimal, x::MOI.VariableIndex
)
    MOI.check_result_index_bounds(model, attr)
    primal_status = MOI.get(model, MOI.PrimalStatus())
    if primal_status == MOI.INFEASIBILITY_CERTIFICATE
        # We claim ownership of the pointer returned by Clp_unboundedRay.
        return _unsafe_wrap_clp_array(
            model.inner,
            Clp_unboundedRay,
            Clp_getNumCols(model.inner),
            x.value;
            own = true,
        )
    elseif primal_status == MOI.FEASIBLE_POINT
        return _unsafe_wrap_clp_array(
            model.inner,
            Clp_getColSolution,
            Clp_getNumCols(model.inner),
            x.value,
        )
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
        # We claim ownership of the pointer returned by Clp_unboundedRay.
        return _unsafe_wrap_clp_array(
            model.inner,
            Clp_unboundedRay,
            Clp_getNumCols(model.inner),
            col_indices;
            own = true,
        )
    elseif primal_status == MOI.FEASIBLE_POINT
        return _unsafe_wrap_clp_array(
            model.inner,
            Clp_getColSolution,
            Clp_getNumCols(model.inner),
            col_indices,
        )
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
    return _unsafe_wrap_clp_array(
        model.inner,
        Clp_getRowActivity,
        Clp_getNumRows(model.inner),
        c.value,
    )
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
    sense = (Clp_getObjSense(model.inner) == -1) ? -1 : 1
    dual_status = MOI.get(model, MOI.DualStatus())
    if dual_status == MOI.FEASIBLE_POINT
        dsol = _unsafe_wrap_clp_array(
            model.inner,
            Clp_getRowPrice,
            Clp_getNumRows(model.inner),
            c.value,
        )
        return sense * dsol
    elseif dual_status == MOI.INFEASIBILITY_CERTIFICATE
        # We claim ownership of the pointer returned by Clp_infeasibilityRay.
        dsol = _unsafe_wrap_clp_array(
            model.inner,
            Clp_infeasibilityRay,
            Clp_getNumRows(model.inner),
            c.value;
            own = true,
        )
        return -sense * dsol
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
    rc = _unsafe_wrap_clp_array(
        model.inner,
        Clp_getReducedCost,
        Clp_getNumCols(model.inner),
        c.value,
    )
    sense = (Clp_getObjSense(model.inner) == -1) ? -1 : 1
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
    rc = _unsafe_wrap_clp_array(
        model.inner,
        Clp_getReducedCost,
        Clp_getNumCols(model.inner),
        c.value,
    )
    sense = (Clp_getObjSense(model.inner) == -1) ? -1 : 1
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
    sense = (Clp_getObjSense(model.inner) == -1) ? -1 : 1
    rc = _unsafe_wrap_clp_array(
        model.inner,
        Clp_getReducedCost,
        Clp_getNumCols(model.inner),
        c.value,
    )
    return sense * rc
end
