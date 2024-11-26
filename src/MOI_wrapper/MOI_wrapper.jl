# Copyright (c) 2013: Clp.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

MOI.Utilities.@product_of_sets(
    _LPProductOfSets,
    MOI.EqualTo{T},
    MOI.LessThan{T},
    MOI.GreaterThan{T},
    MOI.Interval{T},
)

const OptimizerCache = MOI.Utilities.GenericModel{
    Float64,
    MOI.Utilities.ObjectiveContainer{Float64},
    MOI.Utilities.VariablesContainer{Float64},
    MOI.Utilities.MatrixOfConstraints{
        Float64,
        MOI.Utilities.MutableSparseMatrixCSC{
            Float64,
            Cint,
            MOI.Utilities.ZeroBasedIndexing,
        },
        MOI.Utilities.Hyperrectangle{Float64},
        _LPProductOfSets{Float64},
    },
}

Base.show(io::IO, ::Type{OptimizerCache}) = print(io, "Clp.OptimizerCache")

const SCALAR_SETS = Union{
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.EqualTo{Float64},
    MOI.Interval{Float64},
}

"""
    Optimizer()

Create a new Optimizer object.

Set optimizer attributes using `MOI.RawOptimizerAttribute` or
`JuMP.set_optimizer_atttribute`.

For a list of supported parameter names, see `Clp.SUPPORTED_PARAMETERS`.

## Example

```julia
using JuMP, Clp
model = JuMP.Model(Clp.Optimizer)
set_optimizer_attribute(model, "LogLevel", 0)
```
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Ptr{Cvoid}
    solver_options::Ptr{Cvoid}
    options_set::Set{String}
    optimize_called::Bool
    solve_time::Float64
    # Work-around for upstream bug in Clp:
    maximumSeconds::Union{Float64,Nothing}

    function Optimizer()
        model = new(
            Clp_newModel(),
            ClpSolve_new(),
            Set{String}(),
            false,
            0.0,
            nothing,
        )
        finalizer(model) do m
            Clp_deleteModel(m)
            ClpSolve_delete(m.solver_options)
            return
        end
        return model
    end
end

function MOI.default_cache(::Optimizer, ::Type{Float64})
    return MOI.Utilities.UniversalFallback(OptimizerCache())
end

Base.cconvert(::Type{Ptr{Cvoid}}, model::Optimizer) = model
function Base.unsafe_convert(::Type{Ptr{Cvoid}}, model::Optimizer)
    return model.inner::Ptr{Cvoid}
end

# ====================
#   empty functions
# ====================

function MOI.is_empty(model::Optimizer)
    # A problem is empty if it has no variable and no linear constraints
    return (Clp_getNumRows(model) == 0) && (Clp_getNumCols(model) == 0)
end

function MOI.empty!(model::Optimizer)
    # Copy parameters from old model into new model
    old_options = Dict(
        key => MOI.get(model, MOI.RawOptimizerAttribute(key)) for
        key in model.options_set
    )
    empty!(model.options_set)
    Clp_deleteModel(model)
    model.inner = Clp_newModel()
    model.optimize_called = false
    model.solve_time = 0.0
    for (key, value) in old_options
        MOI.set(model, MOI.RawOptimizerAttribute(key), value)
    end
    # Work-around for maximumSeconds
    Clp_setMaximumSeconds(model, something(model.maximumSeconds, -1.0))
    return
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Clp"

function MOI.get(::Optimizer, ::MOI.SolverVersion)
    return "$(Clp_VersionMajor()).$(Clp_VersionMinor()).$(Clp_VersionRelease())"
end

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

function MOI.supports(::Optimizer, param::MOI.RawOptimizerAttribute)
    return param.name in SUPPORTED_PARAMETERS
end

function MOI.set(model::Optimizer, param::MOI.RawOptimizerAttribute, value)
    name = String(param.name)
    push!(model.options_set, name)
    if name == "PrimalTolerance"
        Clp_setPrimalTolerance(model, value)
    elseif name == "DualTolerance"
        Clp_setDualTolerance(model, value)
    elseif name == "DualObjectiveLimit"
        Clp_setDualObjectiveLimit(model, value)
    elseif name == "MaximumIterations"
        Clp_setMaximumIterations(model, value)
    elseif name == "MaximumSeconds"
        Clp_setMaximumSeconds(model, value)
    elseif name == "LogLevel"
        Clp_setLogLevel(model, value)
    elseif name == "Scaling"
        Clp_scaling(model, value)
    elseif name == "Perturbation"
        Clp_setPerturbation(model, value)
    elseif name == "Algorithm"
        Clp_setAlgorithm(model, value)
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

function MOI.get(model::Optimizer, param::MOI.RawOptimizerAttribute)
    name = String(param.name)
    if name == "PrimalTolerance"
        return Clp_primalTolerance(model)
    elseif name == "DualTolerance"
        return Clp_dualTolerance(model)
    elseif name == "DualObjectiveLimit"
        return Clp_dualObjectiveLimit(model)
    elseif name == "MaximumIterations"
        # TODO(odow): Clp doesn't follow convention with maximumIterations!
        return maximumIterations(model)
    elseif name == "MaximumSeconds"
        return Clp_maximumSeconds(model)
    elseif name == "LogLevel"
        return Clp_logLevel(model)
    elseif name == "Scaling"
        return Clp_scalingFlag(model)
    elseif name == "Perturbation"
        return Clp_perturbation(model)
    elseif name == "Algorithm"
        return Clp_algorithm(model)
    elseif name == "PresolveType"
        return ClpSolve_getPresolveType(model.solver_options)
    elseif name == "SolveType"
        return ClpSolve_getSolveType(model.solver_options)
    elseif name == "InfeasibleReturn"
        return ClpSolve_infeasibleReturn(model.solver_options)
    end
    return throw(MOI.UnsupportedAttribute(param))
end

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    push!(model.options_set, "LogLevel")
    Clp_setLogLevel(model, value ? 0 : 1)
    return
end

MOI.get(model::Optimizer, ::MOI.Silent) = Clp_logLevel(model) == 0

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value::Real)
    push!(model.options_set, "MaximumSeconds")
    float_value = convert(Float64, value)
    Clp_setMaximumSeconds(model, float_value)
    model.maximumSeconds = float_value
    return
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    delete!(model.options_set, "MaximumSeconds")
    Clp_setMaximumSeconds(model, -1.0)
    model.maximumSeconds = nothing
    return
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    # TODO(odow): replace with `Clp_maximumSeconds(model)` when upstream
    # is fixed.
    return model.maximumSeconds
end

MOI.supports(::Optimizer, ::MOI.NumberOfThreads) = false

# ========================================
#   Supported constraints and objectives
# ========================================

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:Union{MOI.VariableIndex,MOI.ScalarAffineFunction{Float64}}},
    ::Type{<:SCALAR_SETS},
)
    return true
end

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
)
    return true
end

# =======================
#   `copy_to` function
# =======================

function _index_map(
    src::OptimizerCache,
    index_map,
    ::Type{F},
    ::Type{S},
) where {F,S}
    inner = index_map.con_map[F, S]
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        row = MOI.Utilities.rows(src.constraints, ci)
        inner[ci] = MOI.ConstraintIndex{F,S}(row)
    end
    return
end

function _index_map(
    src::OptimizerCache,
    index_map,
    F::Type{MOI.VariableIndex},
    ::Type{S},
) where {S}
    inner = index_map.con_map[F, S]
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        col = index_map[MOI.VariableIndex(ci.value)].value
        inner[ci] = MOI.ConstraintIndex{F,S}(col)
    end
    return
end

"""
    _index_map(src::OptimizerCache)

Create an `IndexMap` mapping the variables and constraints in `OptimizerCache`
to their corresponding 1-based columns and rows.
"""
function _index_map(src::OptimizerCache)
    index_map = MOI.IndexMap()
    for (i, x) in enumerate(MOI.get(src, MOI.ListOfVariableIndices()))
        index_map[x] = MOI.VariableIndex(i)
    end
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        _index_map(src, index_map, F, S)
    end
    return index_map
end

function MOI.copy_to(dest::Optimizer, src::OptimizerCache)
    @assert MOI.is_empty(dest)
    A = src.constraints.coefficients
    row_bounds = src.constraints.constants
    obj =
        MOI.get(src, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    c = zeros(A.n)
    for term in obj.terms
        c[term.variable.value] += term.coefficient
    end
    Clp_setObjectiveOffset(dest, -obj.constant)
    Clp_loadProblem(
        dest,
        A.n,
        A.m,
        A.colptr,
        A.rowval,
        A.nzval,
        src.variables.lower,
        src.variables.upper,
        c,
        row_bounds.lower,
        row_bounds.upper,
    )
    sense = MOI.get(src, MOI.ObjectiveSense())
    if sense == MOI.MIN_SENSE
        Clp_setObjSense(dest, 1)
    elseif sense == MOI.MAX_SENSE
        Clp_setObjSense(dest, -1)
    else
        @assert sense == MOI.FEASIBILITY_SENSE
        Clp_setObjSense(dest, 0)
    end
    return _index_map(src)
end

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.Utilities.UniversalFallback{OptimizerCache},
)
    MOI.Utilities.throw_unsupported(src)
    return MOI.copy_to(dest, src.model)
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    cache = OptimizerCache()
    src_cache = MOI.copy_to(cache, src)
    cache_dest = MOI.copy_to(dest, cache)
    index_map = MOI.IndexMap()
    for (src_x, cache_x) in src_cache.var_map
        index_map[src_x] = cache_dest[cache_x]
    end
    for (src_ci, cache_ci) in src_cache.con_map
        index_map[src_ci] = cache_dest[cache_ci]
    end
    return index_map
end

# ===============================
#   Optimize and post-optimize
# ===============================

function MOI.optimize!(model::Optimizer)
    t = time()
    Clp_initialSolveWithOptions(model, model.solver_options)
    model.solve_time = time() - t
    model.optimize_called = true
    return
end

MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = model.solve_time

MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = Clp_getNumCols(model)

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    if MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
        N = Clp_getNumCols(model)
        x = _unsafe_wrap_clp_array(model, Clp_unboundedRay, N; own = true)
        c = _unsafe_wrap_clp_array(model, Clp_getObjCoefficients, N)
        return c' * x
    end
    return Clp_getObjValue(model)
end

function _active_bound(sense, d, l, u)
    if -1e100 < l && u < 1e100
        return ifelse(sense * d >= 0, l, u)
    elseif l > -1e100
        return l
    elseif u < 1e100
        return u
    else
        return 0.0
    end
end

function MOI.get(model::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    has_dual_ray =
        MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    M = Clp_getNumRows(model)
    dual_objective_value = 0.0
    if !has_dual_ray
        dual_objective_value -= Clp_objectiveOffset(model)
    end
    row_dual = if has_dual_ray
        _unsafe_wrap_clp_array(model, Clp_infeasibilityRay, M; own = true)
    else
        _unsafe_wrap_clp_array(model, Clp_getRowPrice, M)
    end
    l = _unsafe_wrap_clp_array(model, Clp_getRowLower, M)
    u = _unsafe_wrap_clp_array(model, Clp_getRowUpper, M)
    sense = Clp_getObjSense(model)
    for i in 1:M
        if has_dual_ray
            b = _active_bound(sense, row_dual[i], l[i], u[i])
            dual_objective_value += -sense * b * row_dual[i]
        else
            dual_objective_value +=
                _active_bound(sense, row_dual[i], l[i], u[i]) * row_dual[i]
        end
    end
    N = Clp_getNumCols(model)
    l = _unsafe_wrap_clp_array(model, Clp_getColLower, N)
    u = _unsafe_wrap_clp_array(model, Clp_getColUpper, N)
    if has_dual_ray
        nnz = Clp_getNumElements(model)
        vbeg = _unsafe_wrap_clp_array(model, Clp_getVectorStarts, N)
        vlen = _unsafe_wrap_clp_array(model, Clp_getVectorLengths, N)
        vind = _unsafe_wrap_clp_array(model, Clp_getIndices, nnz)
        vval = _unsafe_wrap_clp_array(model, Clp_getElements, nnz)
        for col in 1:N
            π = 0.0
            for i in vbeg[col] .+ (1:vlen[col])
                π += row_dual[vind[i]+1] * vval[i]
            end
            dual_objective_value +=
                sense * _active_bound(sense, π, l[col], u[col]) * π
        end
    else
        π = _unsafe_wrap_clp_array(model, Clp_getReducedCost, N)
        for col in 1:N
            dual_objective_value +=
                _active_bound(sense, π[col], l[col], u[col]) * π[col]
        end
    end
    return dual_objective_value
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if !model.optimize_called
        return MOI.OPTIMIZE_NOT_CALLED
    end
    st = Clp_status(model)
    if st == 0
        return MOI.OPTIMAL
    elseif st == 1
        return MOI.INFEASIBLE
    elseif st == 2
        return MOI.DUAL_INFEASIBLE
    elseif st == 3
        # No more granular information that "some limit is reached"
        return MOI.OTHER_LIMIT
    end
    return MOI.OTHER_ERROR
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    if !model.optimize_called
        return "MOI.OPTIMIZE_NOT_CALLED"
    end
    st = Clp_status(model)
    if st == 0
        return "0 - optimal"
    elseif st == 1
        return "1 - primal infeasible"
    elseif st == 2
        return "2 - dual infeasible"
    elseif st == 3
        return "3 - stopped on iterations etc"
    else
        @assert st == 4
        return "4 - stopped due to errors"
    end
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    if Clp_primalFeasible(model) != 0
        return 1
    elseif Clp_dualFeasible(model) != 0
        return 1
    elseif Clp_isProvenPrimalInfeasible(model) != 0
        return 1
    elseif Clp_isProvenDualInfeasible(model) != 0
        return 1
    end
    return 0
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    elseif Clp_isProvenDualInfeasible(model) != 0
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif Clp_primalFeasible(model) != 0
        return MOI.FEASIBLE_POINT
    end
    return MOI.UNKNOWN_RESULT_STATUS
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    if attr.result_index != 1
        return MOI.NO_SOLUTION
    elseif Clp_isProvenPrimalInfeasible(model) != 0
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif Clp_dualFeasible(model) != 0
        return MOI.FEASIBLE_POINT
    end
    return MOI.UNKNOWN_RESULT_STATUS
end

# ===================
#   Primal solution
# ===================

function _unsafe_wrap_clp_array(
    model::Optimizer,
    f::Function,
    n::Integer,
    indices = nothing;
    own::Bool = false,
)
    p = f(model)
    if p == C_NULL
        return map(x -> NaN, indices)
    end
    x = unsafe_wrap(Array, p, (n,); own = own)
    return indices === nothing ? x : x[indices]
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    x::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(model, attr)
    primal_status = MOI.get(model, MOI.PrimalStatus())
    if primal_status == MOI.INFEASIBILITY_CERTIFICATE
        # We claim ownership of the pointer returned by Clp_unboundedRay.
        return _unsafe_wrap_clp_array(
            model,
            Clp_unboundedRay,
            Clp_getNumCols(model),
            x.value;
            own = true,
        )
    elseif primal_status == MOI.FEASIBLE_POINT
        return _unsafe_wrap_clp_array(
            model,
            Clp_getColSolution,
            Clp_getNumCols(model),
            x.value,
        )
    end
    return error("Primal solution not available")
end

function MOI.get(
    model::Optimizer,
    attr::MOI.VariablePrimal,
    xs::Vector{MOI.VariableIndex},
)
    MOI.check_result_index_bounds(model, attr)
    col_indices = [idx.value for idx in xs]
    primal_status = MOI.get(model, MOI.PrimalStatus())
    if primal_status == MOI.INFEASIBILITY_CERTIFICATE
        # We claim ownership of the pointer returned by Clp_unboundedRay.
        return _unsafe_wrap_clp_array(
            model,
            Clp_unboundedRay,
            Clp_getNumCols(model),
            col_indices;
            own = true,
        )
    elseif primal_status == MOI.FEASIBLE_POINT
        return _unsafe_wrap_clp_array(
            model,
            Clp_getColSolution,
            Clp_getNumCols(model),
            col_indices,
        )
    end
    return error("Primal solution not available")
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:SCALAR_SETS},
)
    MOI.check_result_index_bounds(model, attr)
    if MOI.get(model, MOI.PrimalStatus()) == MOI.INFEASIBILITY_CERTIFICATE
        return MOI.Utilities.get_fallback(model, attr, c)
    end
    return _unsafe_wrap_clp_array(
        model,
        Clp_getRowActivity,
        Clp_getNumRows(model),
        c.value,
    )
end

# TODO: What happens if model is unbounded / infeasible?
function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.VariableIndex,<:SCALAR_SETS},
)
    MOI.check_result_index_bounds(model, attr)
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end

# =================
#   Dual solution
# =================

# If sense is maximize, we negate all the duals to follow MOI conventions
# Feasibility problems are treated as a minimization

"""
    _farkas_variable_dual(model::Optimizer, col::Cint)

Return a Farkas dual associated with the variable bounds of `col`.

Compute the Farkas dual as:

    ā * x = λ' * A * x <= λ' * b = -β + sum(āᵢ * Uᵢ | āᵢ < 0) + sum(āᵢ * Lᵢ | āᵢ > 0)

The Farkas dual of the variable is ā, and it applies to the upper bound if ā < 0,
and it applies to the lower bound if ā > 0.
"""
function _farkas_variable_dual(model::Optimizer, col::Integer)
    m, n = Clp_getNumRows(model), Clp_getNumCols(model)
    nnz = Clp_getNumElements(model)
    vbeg = _unsafe_wrap_clp_array(model, Clp_getVectorStarts, n)
    vlen = _unsafe_wrap_clp_array(model, Clp_getVectorLengths, n)
    indices = vbeg[col] .+ (1:vlen[col])
    vind = _unsafe_wrap_clp_array(model, Clp_getIndices, nnz, indices)
    vval = _unsafe_wrap_clp_array(model, Clp_getElements, nnz, indices)
    # We need to claim ownership of the pointer returned by Clp_infeasibilityRay.
    λ = _unsafe_wrap_clp_array(model, Clp_infeasibilityRay, m; own = true)
    return sum(v * λ[i+1] for (i, v) in zip(vind, vval))
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:SCALAR_SETS},
)
    MOI.check_result_index_bounds(model, attr)
    n = Clp_getNumRows(model)
    dual_status = MOI.get(model, MOI.DualStatus())
    sense = Clp_getObjSense(model)
    if dual_status == MOI.FEASIBLE_POINT
        dsol = _unsafe_wrap_clp_array(model, Clp_getRowPrice, n, c.value)
        return sense * dsol
    elseif dual_status == MOI.INFEASIBILITY_CERTIFICATE
        # We claim ownership of the pointer returned by Clp_infeasibilityRay.
        return -_unsafe_wrap_clp_array(
            model,
            Clp_infeasibilityRay,
            n,
            c.value;
            own = true,
        )
    end
    return error("Dual solution not available")
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}},
)
    MOI.check_result_index_bounds(model, attr)
    if MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
        return min(0.0, _farkas_variable_dual(model, c.value))
    end
    rc = _unsafe_wrap_clp_array(
        model,
        Clp_getReducedCost,
        Clp_getNumCols(model),
        c.value,
    )
    return min(0.0, Clp_getObjSense(model) * rc)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}},
)
    MOI.check_result_index_bounds(model, attr)
    if MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
        return max(0.0, _farkas_variable_dual(model, c.value))
    end
    rc = _unsafe_wrap_clp_array(
        model,
        Clp_getReducedCost,
        Clp_getNumCols(model),
        c.value,
    )
    return max(0.0, Clp_getObjSense(model) * rc)
end

function MOI.get(
    model::Optimizer,
    attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{
        MOI.VariableIndex,
        <:Union{MOI.Interval{Float64},MOI.EqualTo{Float64}},
    },
)
    MOI.check_result_index_bounds(model, attr)
    if MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
        return _farkas_variable_dual(model, c.value)
    end
    rc = _unsafe_wrap_clp_array(
        model,
        Clp_getReducedCost,
        Clp_getNumCols(model),
        c.value,
    )
    return Clp_getObjSense(model) * rc
end

# Corresponds to the `Status` enum defined in
# https://github.com/coin-or/Clp/blob/8419e63/Clp/src/ClpSimplex.hpp#L114
const _CLP_BASIS_STATUS = Dict(
    # isFree
    0x00 => MOI.BASIC,
    # basic
    0x01 => MOI.BASIC,
    # atUpperBound
    0x02 => MOI.NONBASIC_AT_UPPER,
    # atLowerBound
    0x03 => MOI.NONBASIC_AT_LOWER,
    # superBasic
    0x04 => MOI.SUPER_BASIC,
    # isFixed
    0x05 => MOI.NONBASIC,
)

function _nonbasic_status(status, ::Type{<:MOI.LessThan})
    return status == MOI.NONBASIC_AT_LOWER ? MOI.BASIC : MOI.NONBASIC
end

function _nonbasic_status(status, ::Type{<:MOI.GreaterThan})
    return status == MOI.NONBASIC_AT_LOWER ? MOI.NONBASIC : MOI.BASIC
end

_nonbasic_status(::Any, ::Type{<:MOI.EqualTo}) = MOI.NONBASIC

_nonbasic_status(status, ::Type{<:MOI.Interval}) = status

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintBasisStatus,
    c::MOI.ConstraintIndex{F,S},
) where {F,S}
    code = Clp_getRowStatus(model, c.value - 1)
    status = _CLP_BASIS_STATUS[code]
    if status == MOI.NONBASIC_AT_UPPER || status == MOI.NONBASIC_AT_LOWER
        return _nonbasic_status(status, S)
    end
    return status
end

function MOI.get(
    model::Optimizer,
    ::MOI.VariableBasisStatus,
    vi::MOI.VariableIndex,
)
    code = Clp_getColumnStatus(model, vi.value - 1)
    return _CLP_BASIS_STATUS[code]
end
