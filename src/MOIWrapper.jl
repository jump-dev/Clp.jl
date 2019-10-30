import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

using .ClpCInterface

# Helper functions
_bounds(s::MOI.GreaterThan{Float64}) = (s.lower, Inf)
_bounds(s::MOI.LessThan{Float64}) = (-Inf, s.upper)
_bounds(s::MOI.EqualTo{Float64}) = (s.value, s.value)
_bounds(s::MOI.Interval{Float64}) = (s.lower, s.upper)

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
    # Inner Clp object
    inner::Clp.ClpModel

    # Clp solve options
    solver_options::Clp.ClpSolve

    # To handle MOI.OPTIMIZE_NOT_CALLED status
    optimize_called::Bool

    function Optimizer(;kwargs...)
        inner_model = Clp.ClpModel()
        solver_options = Clp.ClpSolve()
        model = new(inner_model, solver_options, false)

        for (key, value) in kwargs
            MOI.set(model, MOI.RawParameter(key), value)
        end

        return model
    end
end

# ====================
#   empty functions
# ====================
function MOI.is_empty(model::Optimizer)
    # A problem is empty if it has no variable and no linear constraint
    return (Clp.get_num_rows(model.inner) == 0) && (Clp.get_num_cols(model.inner) == 0)
end

function MOI.empty!(model::Optimizer)    
    old_model = model.inner

    # Create new Clp object
    model.inner = Clp.ClpModel()

    # Copy parameters from old model into new model
    for (option, (getter, setter)) in CLP_OPTION_MAP
        value = getter(old_model)
        setter(model.inner, value)
    end

    model.optimize_called = false

    # Free old Clp object
    Clp.ClpCInterface.delete_model(old_model)

    return nothing
end


MOI.get(::Optimizer, ::MOI.SolverName) = "Clp"

# TODO: improve type-stability of `MOI.RawParameter`-related methods.
MOI.supports(::Optimizer, param::MOI.RawParameter) = true

function MOI.set(model::Optimizer, param::MOI.RawParameter, value)
    
    if haskey(CLP_OPTION_MAP, param.name)
        CLP_OPTION_MAP[param.name][2](model.inner, value)
    elseif haskey(SOLVE_OPTION_MAP, param.name)
        SOLVE_OPTION_MAP[param.name][2](model.solver_options, value)
    else
        throw(MOI.UnsupportedAttribute(param))
    end

    return nothing
end

function MOI.get(model::Optimizer, param::MOI.RawParameter)
    if haskey(CLP_OPTION_MAP, param.name)
        return CLP_OPTION_MAP[param.name][1](model.inner)
    elseif haskey(SOLVE_OPTION_MAP, param.name)
        return SOLVE_OPTION_MAP[param.name][1](model.solver_options)
    else
        throw(MOI.UnsupportedAttribute(param))
    end
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    if value
        Clp.set_log_level(model.inner, 0)
    else
        Clp.set_log_level(model.inner, 1)
    end
    return nothing
end

function MOI.get(model::Optimizer, ::MOI.Silent)
    return Clp.log_level(model.inner) == 0
end

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value)
    if value === nothing
        # De-activate time limit
        Clp.set_maximum_seconds(model.inner, -1.0)
    else
        Clp.set_maximum_seconds(model.inner, value)
    end
    return nothing
end

# Clp behaves weird here
MOI.get(model::Optimizer, ::MOI.TimeLimitSec) = Clp.maximum_seconds(model.inner)

MOI.supports(::Optimizer, ::MOI.NumberOfThreads) = false


# ========================================
#   Supported constraints and objectives
# ========================================
MOI.supports_constraint(::Optimizer,
    ::Type{MOI.ScalarAffineFunction{Float64}}, ::Type{<:SCALAR_SETS}
) = true

MOI.supports_constraint(::Optimizer,
    ::Type{MOI.SingleVariable}, ::Type{<:SCALAR_SETS}
) = true

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

MOI.supports(::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}
) = true

# =======================
#   `copy_to` function
# =======================

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.ModelLike;
    copy_names=false
)
    # Check that all constraints and objectives
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        MOI.supports_constraint(dest, F, S) || throw(MOI.UnsupportedConstraint{F, S}(
            "Clp.Optimizer does not support constraints of type $F-in-$S."
        ))
    end
    fobj_type = MOI.get(src, MOI.ObjectiveFunctionType())
    MOI.supports(dest, MOI.ObjectiveFunction{fobj_type}()) || throw(
        MOI.UnsupportedAttribute(MOI.ObjectiveFunction(fobj))
    )

    # Empty dest
    MOI.empty!(dest)

    # Maps the indices of the 
    mapping = MOIU.IndexMap()

    # First create variables (including bounds)

    # Create variables
    # TODO: it may be more efficient to create the problem all at once
    x_src = MOI.get(src, MOI.ListOfVariableIndices())
    for (j, x) in enumerate(x_src)
        # Variable j corresponds to j-th variable of x_src
        mapping.varmap[x] = MOI.VariableIndex(j)
        
        # Possible bound combinations:
        #=
            * No bound (default case)
            * Interval
            * EqualTo
            * LessThan
            * GreaterThan
            * LessThan && GreaterThan
        =#
        lb, ub = -Inf, Inf  # Default case: free variable

        idx_RG = MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}(x.value)
        idx_ET = MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}(x.value)
        idx_LT = MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}(x.value)
        idx_GT = MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}(x.value)
        if MOI.is_valid(src, idx_LT)
            # Get the (upper) bound
            s = MOI.get(src, MOI.ConstraintSet(), idx_LT)
            ub = s.upper

            # Track index of upper-bound constraint
            mapping.conmap[idx_LT] = MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}(j)
        end
        if MOI.is_valid(src, idx_GT)
            # Get the (lower) bound
            s = MOI.get(src, MOI.ConstraintSet(), idx_GT)
            lb = s.lower

            # Track index of upper-bound constraint
            mapping.conmap[idx_GT] = MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}(j)
        end
        if MOI.is_valid(src, idx_RG)
            # Get upper and lower bound
            s = MOI.get(src, MOI.ConstraintSet(), idx_RG)
            lb, ub = s.lower, s.upper

            # Track index of upper-bound constraint
            mapping.conmap[idx_RG] = MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}(j)
        end
        if MOI.is_valid(src, idx_ET)
            # Get the bounds
            s = MOI.get(src, MOI.ConstraintSet(), idx_ET)
            lb, ub = s.value, s.value

            # Track index of upper-bound constraint
            mapping.conmap[idx_ET] = MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Float64}}(j)
        end

        Clp.add_column(dest.inner, Cint(0), Int32[], Float64[], lb, ub, 0.0)
    end

    # Create linear constraints
    nrows = 0  # number of rows in the problem
    for (F, S) in MOI.get(src, MOI.ListOfConstraints())
        # Check if constraint is a linear constraint
        (F == MOI.ScalarAffineFunction{Float64}) || continue

        for index in MOI.get(src, MOI.ListOfConstraintIndices{F, S}())
            # Retrieve problem data
            # TODO: ensure that `(f, s)` are in canonical form
            f = MOIU.canonical(MOI.get(src, MOI.ConstraintFunction(), index))
            s = MOI.get(src, MOI.ConstraintSet(), index)

            # Extract bounds
            # Constraint writes a'x + a0 - in - set
            # So we offset the bounds by -a0
            lb, ub = _bounds(s) .- f.constant

            # Identify constraint terms
            coeffs = [t.coefficient for t in f.terms]
            cols::Vector{Int} = [
                mapping.varmap[t.variable_index].value
                for t in f.terms
            ]

            # Add constraint
            Clp.add_row(dest.inner, Cint(length(f.terms)), Cint.(cols .- 1), coeffs, lb, ub)

            # Track index
            mapping.conmap[index] = MOI.ConstraintIndex{F, S}(nrows + 1)
            nrows += 1
        end
    end

    # Objective function
    sense = MOI.get(src, MOI.ObjectiveSense())
    if sense == MOI.MIN_SENSE
        Clp.set_obj_sense(dest.inner, 1)   # minimize
    elseif sense == MOI.MAX_SENSE
        Clp.set_obj_sense(dest.inner, -1)  # maximize
    elseif sense == MOI.FEASIBILITY_SENSE
        Clp.set_obj_sense(dest.inner, 0)   # feasibility
    else
        # Should not be reached, but just in case
        error("Unexpected optimization sense $sense")
    end

    fobj = MOI.get(src, MOI.ObjectiveFunction{fobj_type}())
    obj_coeffs = zeros(length(x_src))
    obj_offset = 0.0

    if fobj_type == MOI.ScalarAffineFunction{Float64}
        # Record objective coeffs
        for term in fobj.terms
            x = mapping.varmap[term.variable_index]
            # Here we do += instead of = to accumulate dupplicates
            obj_coeffs[x.value] += term.coefficient
        end
        obj_offset = fobj.constant
    else
        # Should not be reached here, but just in case
        error("Unexpected objective $fobj")
    end

    # Set objective coefficients
    Clp.chg_obj_coefficients(dest.inner, obj_coeffs)
    # Set objective offset
    # Clp seems to negates the objective offset
    Clp.set_objective_offset(dest.inner, -obj_offset)
    
    return mapping
end

# ===============================
#   Optimize and post-optimize
# ===============================

function MOI.optimize!(model::Optimizer)
    Clp.initial_solve_with_options(model.inner, model.solver_options)
    model.optimize_called = true
    return nothing
end

MOI.get(model::Optimizer, ::MOI.NumberOfVariables) = Clp.get_num_cols(model.inner)

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    # Check result index
    MOI.check_result_index_bounds(model, attr)

    return Clp.get_obj_value(model.inner)
end

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)

    model.optimize_called || return MOI.OPTIMIZE_NOT_CALLED

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
    model.optimize_called || return "MOI.OPTIMIZE_NOT_CALLED"

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
    st = MOI.get(model, MOI.TerminationStatus())

    # TODO: check that this is correct
    if st == MOI.OPTIMIZE_NOT_CALLED || st == MOI.OTHER_ERROR
        return 0
    end
    return 1
end

# Solution status
function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    # Check result index
    attr.N == 1 || return MOI.NO_SOLUTION

    if Clp.is_proven_dual_infeasible(model.inner)
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif Clp.primal_feasible(model.inner)
        return MOI.FEASIBLE_POINT
    else
        return MOI.UNKNOWN_RESULT_STATUS
    end
end

function MOI.get(model::Optimizer, attr::MOI.DualStatus)
    # Check result index
    attr.N == 1 || return MOI.NO_SOLUTION

    if Clp.is_proven_primal_infeasible(model.inner)
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
function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, x::MOI.VariableIndex)
    # Check result index
    MOI.check_result_index_bounds(model, attr)

    # We assume the variable index is valid,
    # since Clp should be accessed via a CachingOptimizer
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

function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, xs::Vector{MOI.VariableIndex})
    # Check result index
    MOI.check_result_index_bounds(model, attr)

    col_indices = [idx.value for idx in xs]
    
    # We assume all variable indices is valid,
    # since Clp should be accessed via a CachingOptimizer
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

function MOI.get(model::Optimizer, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}
) where{S <:SCALAR_SETS}
    # Check result index
    MOI.check_result_index_bounds(model, attr)

    # We assume that the constraint index is valid,
    # since Clp should be accessed via a CachingOptimizer
    # TODO: What happens if model is unbounded / infeasible?

    # Get primal row solution
    sol = Clp.get_row_activity(model.inner)
    return sol[c.value]
end

function MOI.get(model::Optimizer, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{S <: SCALAR_SETS}
    # Check result index
    MOI.check_result_index_bounds(model, attr)

    # TODO: What happens if model is unbounded / infeasible?
    return MOI.get(model, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end


# =================
#   Dual solution
# =================
# If sense is maximize, we negate all the duals to follow MOI conventions
# Feasibility problems are treated as a minimization
function MOI.get(model::Optimizer, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, S}
) where{S <: SCALAR_SETS}
    # Check result index
    MOI.check_result_index_bounds(model, attr)

    sense = (Clp.get_obj_sense(model.inner) == -1) ? -1 : 1

    dual_status = MOI.get(model, MOI.DualStatus())
    if dual_status == MOI.FEASIBLE_POINT
        # Return dual solution
        dsol = Clp.get_row_price(model.inner)
        return sense * dsol[c.value]
    elseif dual_status == MOI.INFEASIBILITY_CERTIFICATE
        # Return infeasibility ray
        dsol = Clp.infeasibility_ray(model.inner)
        return sense * dsol[c.value]
    else
        # Throw error if dual solution is not available
        error("Dual solution not available")
    end
end

# Reduced costs
# TODO: what happens if problem is unbounded / infeasible?
function MOI.get(model::Optimizer, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Float64}}   
)
    # Check result index
    MOI.check_result_index_bounds(model, attr)

    rc = Clp.get_reduced_cost(model.inner)[c.value]
    sense = (Clp.get_obj_sense(model.inner) == -1) ? -1 : 1

    # Dual should be non-positive
    if sense == 1 && rc <= 0.0
        return rc
    elseif sense == -1 && rc >= 0.0
        return -rc
    else
        return 0.0
    end
end

function MOI.get(model::Optimizer, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Float64}}   
)
    # Check result index
    MOI.check_result_index_bounds(model, attr)

    rc = Clp.get_reduced_cost(model.inner)[c.value]
    sense = (Clp.get_obj_sense(model.inner) == -1) ? -1 : 1

    # Dual should be non-negative   
    if sense == 1 && rc >= 0.0
        return rc
    elseif sense == -1 && rc <= 0.0
        return -rc
    else
        return 0.0
    end
end

function MOI.get(model::Optimizer, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}   
) where{S <: Union{MOI.Interval{Float64}, MOI.EqualTo{Float64}}}
    # Check result index
    MOI.check_result_index_bounds(model, attr)
    
    sense = (Clp.get_obj_sense(model.inner) == -1) ? -1 : 1
    return sense * Clp.get_reduced_cost(model.inner)[c.value]
end