export ClpOptimizer

using LinQuadOptInterface
using Clp.ClpCInterface

function replaceInf(x)
    for i in 1:length(x)
        if x[i] > 1e20
            x[i] = Inf
        elseif x[i] < -1e20
            x[i] = -Inf
        end
    end
    return x
end

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
    (LQOI.SinVar, MOI.ZeroOne),
    (LQOI.SinVar, MOI.Integer),
    (LQOI.VecVar, MOI.Nonnegatives),
    (LQOI.VecVar, MOI.Nonpositives),
    (LQOI.VecVar, MOI.Zeros),
    (LQOI.VecLin, MOI.Nonnegatives),
    (LQOI.VecLin, MOI.Nonpositives),
    (LQOI.VecLin, MOI.Zeros)
]

mutable struct ClpOptimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    # env
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
   #:Algorithm => set_algorithm
   )
# These options are set by using the ClpSolve object
const solveoptionmap = Dict(
   :PresolveType => set_presolve_type,
   :SolveType => set_solve_type,
   :InfeasibleReturn => set_infeasible_return,
   )

function setoption(m::ClpOptimizer, name::Symbol, value)
    if haskey(optionmap, name)
        optionmap[name](m.inner,value)
    elseif haskey(solveoptionmap, name)
        solveoptionmap[name](m.solveroptions,value)
    else
        error("Unrecognized option: $name")
    end
end

function ClpOptimizer(;kwargs...)
    # env = Env()
    m = ClpOptimizer(nothing)
    # m.env = env
    m.params = Dict{String,Any}()
    MOI.empty!(m)
    for (name,value) in kwargs
        m.params[Symbol(name)] = value
    end
    return m
end

LQOI.LinearQuadraticModel(::Type{ClpOptimizer},env) = ClpModel() 

LQOI.supported_constraints(s::ClpOptimizer) = SUPPORTED_CONSTRAINTS
LQOI.supported_objectives(s::ClpOptimizer) = SUPPORTED_OBJECTIVES

"""
Change the bounds of the variable. The sense of the upperbound
is given by `backend_type(m, Val{:Upperbound}())`. The sense
of the lowerbound is given by `backend_type(m, Val{:Lowerbound}())`
"""
function LQOI.change_variable_bounds!(instance::ClpOptimizer, cols::Vector{Int}, values::Vector{Float64}, 
        senses::Vector)
    
    upperbounds = replaceInf(get_col_upper(instance.inner))    
    lowerbounds = replaceInf(get_col_lower(instance.inner))    
    
    for i in 1:length(senses)
        if senses[i] == Cchar('U')
            upperbounds[cols[i]] = values[i]
        elseif senses[i] == Cchar('L')
            lowerbounds[cols[i]] = values[i]
        else            
            error("sense is " * string(senses[i]) * ", but only " * Cchar('U') *  
            " and " * Cchar('L') * " senses are supported")
        end
    end
    
    chg_column_upper(instance.inner, upperbounds)
    chg_column_lower(instance.inner, lowerbounds)
end

"""
get_variable_lowerbound(m, col::Int)::Float64

Get the lower bound of the variable in 1-indexed column `col` of the model `m`.
"""
function LQOI.get_variable_lowerbound(instance::ClpOptimizer, col::Int)::Float64
    lower = replaceInf(get_col_lower(instance.inner))
    return lower[col];  
end

"""
get_variable_upperbound(m, col::Int)::Float64

Get the upper bound of the variable in 1-indexed column `col` of the model `m`.
"""
function LQOI.get_variable_upperbound(instance::ClpOptimizer, col::Int)::Float64
    upper = replaceInf(get_col_upper(instance.inner))
    return upper[col];  
end

"""
get_number_linear_constraints(m)::Int

Get the number of linear constraints in the model `m`.
"""
function LQOI.get_number_linear_constraints(instance::ClpOptimizer)
    return get_num_rows(instance.inner)
end

"""
    add_linear_constraints!(m, A::CSRMatrix{Float64}, sense::Vector{Cchar}, rhs::Vector{Float64})::Void

Adds linear constraints of the form `Ax (sense) rhs` to the model `m`.
`sense` and `rhs` contain one element for each row in `A`.
The `sense` is given by `backend_type(m, set)`.
Ranged constraints (`set=MOI.Interval`) should be added via `add_ranged_constraint!` instead.
See also: `LinQuadOptInterface.CSRMatrix`.
"""
# function LQOI.add_linear_constraints!(instance::ClpOptimizer, rows::Vector{Int}, cols::Vector{Int},
#         coefs::Vector{Float64}, sense::Vector{Cchar}, rhs::Vector{Float64})::Void
function LQOI.add_linear_constraints!(instance::ClpOptimizer, A::LQOI.CSRMatrix{Float64}, sense::Vector{Cchar}, 
        rhs::Vector{Float64})::Void
    
    rows = A.row_pointers
    cols = A.columns
    coefs = A.coefficients
           
    nbrows = length(rhs)
    for row in 1:nbrows
        elements = Vector{Float64}()
        columns = Vector{Int32}()
        lower = -Inf
        upper = Inf
        
        if sense[row] == Cchar('L')
            upper = rhs[row]
        elseif sense[row] == Cchar('G')
            lower = rhs[row]
        elseif sense[row] == Cchar('E')
            upper = lower = rhs[row]
        else
            error("sense must be Cchar(x) where x is in ['L','G',E']")
        end
        
        number_in_row = 0
        first = rows[row]
        last = (row==length(sense)) ? length(cols) : rows[row+1]-1
        for i in first:last
            number_in_row += 1
            push!(elements, coefs[i])
            push!(columns, cols[i] - 1)
        end  
        
        add_row(instance.inner, Cint(number_in_row), columns, elements, lower, upper)
    end
end

# 
# """
# add_ranged_constraint!(m, rows::Vector{Int}, cols::Vector{Int},
#     coefs::Vector{Float64}, lowerbound::Vector{Float64}, upperbound::Vector{Float64})
# 
# Adds linear constraints of the form `lowerbound <= Ax <= upperbound` to the
# model `m`.
# 
# The A matrix is given in triplet form `A[rows[i], cols[i]] = coef[i]` for all
# `i`,
# 
# This is a special case compared to standard `add_linear_constraints!` since it
# is often implemented via multiple API calls.
# """
# # function add_ranged_constraints! end
# 
# """
# modify_ranged_constraint!(m, rows::Vector{Int}, lowerbound::Vector{Float64}, upperbound::Vector{Float64})
# 
# Modify the lower and upperbounds of a ranged constraint in the model `m`.
# 
# This is a special case compared to standard the `change_coefficient!` since it
# is often implemented via multiple API calls.
# """
# # function modify_ranged_constraints! end
# 
# 

"""
get_rhs(m, row::Int)::Float64

Get the right-hand side of the linear constraint in the 1-indexed row `row` in
the model `m`.
"""
function LQOI.get_rhs(instance::ClpOptimizer, row::Int)::Float64
    lower = replaceInf(get_row_lower(instance.inner))
    upper = replaceInf(get_row_upper(instance.inner))
    lb = lower[row]
    ub = upper[row]
    
    if lb > -Inf
        return lb
    elseif ub < Inf
        return ub
    else
        error("Either row_lower or row_upper must be of abs less than 1e20")
    end
end

"""
get_linear_constraint(m, row::Int)::Tuple{Vector{Int}, Vector{Float64}}

Get the linear component of the constraint in the 1-indexed row `row` in
the model `m`. Returns a tuple of `(cols, vals)`.
"""
function LQOI.get_linear_constraint(instance::ClpOptimizer, row::Int)::Tuple{Vector{Int}, Vector{Float64}}
    A = get_constraint_matrix(instance.inner)
    A_row = A[row,:]
    return (Array{Int}(A_row.nzind), A_row.nzval)
end

"""
change_objective_coefficient!(m, col, coef)

Set the linear coefficient of the variable in column `col` to `coef` in the objective function.
"""
function LQOI.change_objective_coefficient!(instance::ClpOptimizer, col, coef)
    objcoefs = get_obj_coefficients(instance.inner)
    objcoefs[col] = coef
    chg_obj_coefficients(instance.inner, objcoefs)
end

"""
change_rhs_coefficient!(m, row, coef)

Set the rhs of the constraint in row `row` to `coef`.
"""
function LQOI.change_rhs_coefficient!(instance::ClpOptimizer, row, coef)
    lower = replaceInf(get_row_lower(instance.inner))
    upper = replaceInf(get_row_upper(instance.inner))
    lb = lower[row]
    ub = upper[row]
    
    if lb > -Inf
        lower[row] = coef
        chg_row_lower(instance.inner, lower)
    elseif ub < Inf
        upper[row] = coef
        chg_row_upper(instance.inner, upper)
    else
        error("Either row_lower or row_upper must be of abs less than 1e20")
    end
end

 
"""
delete_linear_constraints!(m, start_row::Int, end_row::Int)::Void

Delete the linear constraints `start_row`, `start_row+1`, ..., `end_row` from
the model `m`.
"""
function LQOI.delete_linear_constraints!(instance::ClpOptimizer, start_row::Int, end_row::Int)::Void    
    which = [Int32(i-1) for i in start_row:end_row]
    delete_rows(instance.inner, which)
end


"""
change_linear_constraint_sense!(m, rows::Vector{Int}, sense::Vector{Cchar})::Void

Change the sense of the linear constraints in `rows` to `sense`.

`sense` is the output of `backend_type(m, set)`, where `set`
is the corresponding set for the row `rows[i]`.

"""
function LQOI.change_linear_constraint_sense!(instance::ClpOptimizer, rows::Vector{Int}, sense::Vector{Cchar} )::Void
    lower = replaceInf(get_row_lower(instance.inner))
    upper = replaceInf(get_row_upper(instance.inner))
        
    for (i,row) in enumerate(rows)                
        lb = lower[row]
        ub = upper[row]
        if lb > -Inf
            rhs = lb
        elseif ub < Inf
            rhs = ub
        else
            error("Either row_lower or row_upper must be of abs less than 1e20")
        end
        
        sense = sense[i]
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



"""
set_linear_objective!(m, cols::Vector{Int}, coefs::Vector{Float64})::Void

Set the linear component of the objective.
"""
function LQOI.set_linear_objective!(instance::ClpOptimizer, cols::Vector{Int}, coefs::Vector{Float64})::Void
    obj_in = zeros(Float64, get_num_cols(instance.inner))
    for i in 1:length(cols)
        obj_in[cols[i]] = coefs[i]
    end
    chg_obj_coefficients(instance.inner, obj_in) 
end

"""
change_objective_sense!(m, sense::Symbol)::Void

Change the optimization sense of the model `m` to `sense`. `sense` must be
`:min` or `:max`.
"""
function LQOI.change_objective_sense!(instance::ClpOptimizer, sense::Symbol)
    if sense == :min
        set_obj_sense(instance.inner, 1.0)
    elseif sense == :max
        set_obj_sense(instance.inner, -1.0)
    else
        error("sense must be either :min or :max")
    end
end


"""
get_linear_objective!(m, x::Vector{Float64})

Change the linear coefficients of the objective and store
in `x`.
"""
function LQOI.get_linear_objective!(instance::ClpOptimizer, x::Vector{Float64})
    obj = get_obj_coefficients(instance.inner)
    copy!(x, obj)
end

# """
# get_objectivesense(m)::MOI.OptimizationSense
# 
# Get the optimization sense of the model `m`.
# """
# # function get_objectivesense end
#

"""
solve_linear_problem!(m)::Void

Solve a linear program `m`.
"""
function LQOI.solve_linear_problem!(instance::ClpOptimizer)::Void
    solveroptions = ClpSolve()
    model = instance.inner
    for (name, value) in instance.params
        if haskey(optionmap, name)
            optionmap[name](model,value)
        elseif haskey(solveoptionmap, name)
            solveoptionmap[name](solveroptions,value)
        else
            error("Unrecognized option: $name")
        end
    end
    
    initial_solve_with_options(instance.inner, solveroptions)
    return
end

"""
get_variable_primal_solution!(m, x::Vector{Float64})

Get the primal solution for the variables in the model `m`, and
store in `x`. `x`must have one element for each variable.
"""
function LQOI.get_variable_primal_solution!(instance::ClpOptimizer, x::Vector{Float64})
    solution = primal_column_solution(instance.inner)
    for i in 1:length(solution)
        x[i] = solution[i]
    end
end

"""
get_linear_primal_solution!(m, x::Vector{Float64})

Given a set of linear constraints `l <= a'x <= b` in the model `m`, get the
constraint primal `a'x` for each constraint, and store in `x`.
`x` must have one element for each linear constraint.
"""
function LQOI.get_linear_primal_solution!(instance::ClpOptimizer, x::Vector{Float64})
    solution = primal_row_solution(instance.inner)
    for i in 1:length(solution)
        x[i] = solution[i]
    end
end

"""
get_variable_dual_solution!(m, x::Vector{Float64})

Get the dual solution (reduced-costs) for the variables in the model `m`, and
store in `x`. `x`must have one element for each variable.
"""
function LQOI.get_variable_dual_solution!(instance::ClpOptimizer, x::Vector{Float64})
    solution = dual_column_solution(instance.inner)
    for i in 1:length(solution)
        x[i] = solution[i]
    end
end

"""
get_linear_dual_solution!(m, x::Vector{Float64})

Get the dual solution for the linear constraints in the model `m`, and
store in `x`. `x`must have one element for each linear constraint.
"""
function LQOI.get_linear_dual_solution!(instance::ClpOptimizer, x::Vector{Float64})
    solution = dual_row_solution(instance.inner)
    for i in 1:length(solution)
        x[i] = solution[i]
    end
end

"""
get_objective_value(m)

Get the objective value of the solved model `m`.
"""
function LQOI.get_objective_value(instance::ClpOptimizer)
    return objective_value(instance.inner)
end

# 
# """
# get_objective_bound(m)
# 
# Get the objective bound of the model `m`.
# """
# # function get_objective_bound end
# 
# """
# get_relative_mip_gap(m)
# 
# Get the relative MIP gap of the solved model `m`.
# """
# # function get_relative_mip_gap end
# 
# """
# get_iteration_count(m)
# 
# Get the number of simplex iterations performed during the most recent
# optimization of the model `m`.
# """
# # function get_iteration_count end
# 
# """
# get_barrier_iterations(m)
# 
# Get the number of barrier iterations performed during the most recent
# optimization of the model `m`.
# """
# # function get_barrier_iterations end
# 
# """
# get_node_count(m)
# 
# Get the number of branch-and-cut nodes expolored during the most recent
# optimization of the model `m`.
# """
# # function get_node_count end
# 
# """
# get_farkas_dual!(m, x::Vector{Float64})
# 
# Get the farkas dual (certificate of primal infeasiblility) for the linear
# constraints in the model `m`, and store in `x`. `x`must have one element for
# each linear constraint.
# """
# # function get_farkas_dual! end
# 
# """
# get_unbounded_ray!(m, x::Vector{Float64})
# 
# Get the unbounded ray (certificate of dual infeasiblility) for the linear
# constraints in the model `m`, and store in `x`. `x`must have one element for
# each variable.
# """
# # function get_unbounded_ray! end
# 

"""
get_termination_status(m)

Get the termination status of the model `m`.
"""
function LQOI.get_termination_status(instance::ClpOptimizer)
    s = ClpCInterface.status(instance.inner)
    if s == 0 
        return MOI.Success
    elseif s == 1
        return MOI.InfeasibleNoResult
    elseif s == 2
        return MOI.UnboundedNoResult
    # if s in [0, 1, 2] 
    #     return MOI.Success
    elseif s == 3
        return MOI.OtherLimit
    elseif s == 4
        return MOI.OtherError
    else
        error("status returned by Clp must be in [0,1,2,3,4]")
    end
end

"""
get_primal_status(m)

Get the primal status of the model `m`.
"""
function LQOI.get_primal_status(instance::ClpOptimizer)
    if primal_feasible(instance.inner)
        return MOI.FeasiblePoint
    else    
        return MOI.UnknownResultStatus
    end
end

"""
get_dual_status(m)

Get the dual status of the model `m`.
"""
function LQOI.get_dual_status(instance::ClpOptimizer)
    if dual_feasible(instance.inner)
        return MOI.FeasiblePoint     
    else
        return MOI.UnknownResultStatus
    end    
end

"""
get_number_variables(m)::Int

Get the number of variables in the model `m`.
"""
function LQOI.get_number_variables(instance::ClpOptimizer)
    return get_num_cols(instance.inner)
end

"""
add_variables!(m, n::Int)::Void

Add `n` new variables to the model `m`.
"""
function LQOI.add_variables!(instance::ClpOptimizer, n::Int)::Void
    for i in 1:n
        numberInColumn = 0
        rows = Vector{Int32}()
        elements = Vector{Float64}()
        add_column(instance.inner, numberInColumn,  rows, elements, -Inf, +Inf, 0.0)
    end
end

"""
delete_variables!(m, start_col::Int, end_col::Int)::Void

Delete the columns `start_col`, `start_col+1`, ..., `end_col` from the model `m`.
"""
function LQOI.delete_variables!(instance::ClpOptimizer, start_col::Int, end_col::Int)::Void
    which = [Int32(i-1) for i in start_col:end_col]
    delete_columns(instance.inner, which)
end