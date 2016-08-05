module ClpCInterface

import Compat: String, unsafe_wrap

# Load binary dependencies via Cbc package
import Cbc

export
    # Types
    CoinBigIndex,
    CoinBigDouble,
    ClpModel,
    ClpSolve,

    # Methods
    load_problem,
    load_quadratic_objective,
    read_mps,
    copy_in_integer_information,
    delete_integer_information,
    resize,
    delete_rows,
    add_rows,
    delete_columns,
    add_columns,
    chg_row_lower,
    chg_row_upper,
    chg_column_lower,
    chg_column_upper,
    chg_obj_coefficients,
    drop_names,
    copy_names,
    get_num_rows,
    get_num_cols,
    number_rows,
    number_cols,
    primal_tolerance,
    set_primal_tolerance,
    dual_tolerance,
    set_dual_tolerance,
    dual_objective_limit,
    set_dual_objective_limit,
    objective_offset,
    set_objective_offset,
    problem_name,
    set_problem_name,
    number_iterations,
    get_iteration_count,
    set_number_iterations,
    maximum_iterations,
    set_maximum_iterations,
    maximum_seconds,
    set_maximum_seconds,
    hit_maximum_iterations,
    status,
    set_problem_status,
    secondary_status,
    set_secondary_status,
    optimization_direction,
    set_optimization_direction,
    get_row_activity,
    get_col_solution,
    primal_column_solution,
    set_col_solution,
    primal_row_solution,
    get_row_price,
    dual_row_solution,
    get_reduced_cost,
    dual_column_solution,
    get_row_lower,
    row_lower,
    get_row_upper,
    row_upper,
    get_obj_coefficients,
    objective,
    get_col_lower,
    column_lower,
    get_col_upper,
    column_upper,
    get_num_elements,
    get_vector_starts,
    get_indices,
    get_constraint_matrix,
    get_vector_lengths,
    get_elements,
    get_obj_value,
    objective_value,
    integer_information,
    infeasibility_ray,
    unbounded_ray,
    statusExists,
    status_array,
    copyin_status,
    get_column_status,
    get_row_status,
    set_column_status,
    set_row_status,
    register_call_back,
    clear_call_back,
    log_level,
    set_log_level,
    length_names,
    row_name,
    column_name,
    initial_solve,
    initial_solve_with_options,
    initial_dual_solve,
    initial_primal_solve,
    initial_barrier_solve,
    initial_barrier_no_cross_solve,
    dual,
    primal,
    idiot,
    scaling,
    scaling_flag,
    crash,
    primal_feasible,
    dual_feasible,
    dual_bound,
    set_dual_bound,
    infeasibility_cost,
    set_infeasibility_cost,
    perturbation,
    set_perturbation,
    algorithm,
    set_algorithm,
    sum_dual_infeasibilities,
    number_dual_infeasibilities,
    sum_primal_infeasibilities,
    number_primal_infeasibilities,
    save_model,
    restore_model,
    check_solution,
    is_abandoned,
    is_proven_optimal,
    is_proven_primal_infeasible,
    is_proven_dual_infeasible,
    is_primal_objective_limit_reached,
    is_dual_objective_limit_reached,
    is_iteration_limit_reached,
    get_obj_sense,
    set_obj_sense,
    print_model,
    get_small_element_value,
    set_small_element_value,

    # ClpSolve methods
    set_special_option,
    get_special_option,
    set_presolve_type,
    get_presolve_type,
    set_solve_type,
    get_solve_type,
    get_presolve_passes,
    get_extra_info,
    set_infeasible_return,
    infeasible_return,
    do_dual,
    set_do_dual,
    do_singleton,
    set_do_singleton,
    do_doubleton,
    set_do_doubleton,
    do_tripleton,
    set_do_tripleton,
    do_tighten,
    set_do_tighten,
    do_forcing,
    set_do_forcing,
    do_implied_free,
    set_do_implied_free,
    do_dupcol,
    set_do_dupcol,
    do_duprow,
    set_do_duprow,
    do_singleton_column,
    set_do_singleton_column,
    presolve_actions,
    set_presolve_actions,
    substitution,
    set_substitution

import Base.pointer


## Shared library interface setup
#{{{

macro clp_ccall(func, args...)
    f = "Clp_$(func)"
    quote
        ccall(($f,Cbc.libclp), $(args...))
    end
end

macro clpsolve_ccall(func, args...)
    f = "ClpSolve_$(func)"
    quote
        ccall(($f,Cbc.libclp), $(args...))
    end
end


# Note: we assume COIN_BIG_INDEX and COIN_BIG_DOUBLE
# were not defined when compiling Clp (which is the
# default)
typealias CoinBigIndex Int32
typealias CoinBigDouble Float64

#}}}

## Main types definitions
#{{{
type ClpModel
    p::Ptr{Void}
    function ClpModel()
        p = @clp_ccall newModel Ptr{Void} ()
        prob = new(p)
        finalizer(prob, delete_model)
        return prob
    end
end

function delete_model(model::ClpModel)
    if model.p == C_NULL
        return
    end
    @clp_ccall deleteModel Void (Ptr{Void},) model.p
    model.p = C_NULL
    return
end

type ClpSolve
    p::Ptr{Void}
    function ClpSolve()
        p = @clpsolve_ccall new Ptr{Void} ()
        prob = new(p)
        finalizer(prob, delete_solve)
        return prob
    end
end

function delete_solve(solve::ClpSolve)
    if solve.p == C_NULL
        return
    end
    @clpsolve_ccall delete Void (Ptr{Void},) solve.p
    solve.p = C_NULL
    return
end

pointer(model::ClpModel) = model.p
#}}}

## Check functions for internal use
#{{{
# Functions which perform all sorts of
# sanity checks on input parameters and
# throw exceptions in case of errors.
# Ideally, it should never be possible
# to pass an invalid parameter to the
# underlying Clp API.

function _jl__check_model(model::ClpModel)
    if model.p == C_NULL
        error("Invalid ClpModel")
    end
    return true
end

function _jl__check_solve(solve::ClpSolve)
    if solve.p == C_NULL
        error("Invalid ClpSolve")
    end
    return true
end

function _jl__check_row_is_valid(model::ClpModel, row::Integer)
    num_rows = @clp_ccall getNumRows Int32 (Ptr{Void},) model.p
    if !(1 <= row <= num_rows)
        error("Invalid row $row (must be 1 <= row <= $num_rows)")
    end
    return true
end

function _jl__check_col_is_valid(model::ClpModel, col::Integer)
    num_cols = @clp_ccall getNumCols Int32 (Ptr{Void},) model.p
    if !(1 <= col <= num_cols)
        error("Invalid col $col (must be 1 <= col <= $num_cols)")
    end
    return true
end

function _jl__check_file_is_readable(filename::String)
    try
        f = open(filename, "r")
        close(f)
    catch err
        error("file $filename not readable")
    end
    return true
end
#}}}

## CLP functions
#{{{

# inspired by GLPK interface
typealias VecOrNothing Union{Vector,Void}
function vec_or_null{T}(::Type{T}, a::VecOrNothing, len::Integer)
    if isequal(a, nothing) || isa(a, Array{Void}) # on 0.3, [] is Array{Void}
        return C_NULL
    else # todo: helpful message if convert fails
        if length(a) != len
            error("Expected vector to have length $len")
        end
        return convert(Vector{T},a)
    end
end

# Load model - loads some stuff and initializes others
# Loads a problem (the constraints on the
# rows are given by lower and upper bounds). If a pointer is NULL then the
# following values are the default:
#
# col_ub: all columns have upper bound infinity
# col_lb: all columns have lower bound 0
# row_ub: all rows have upper bound infinity
# row_lb: all rows have lower bound -infinity
# obj: all variables have 0 objective coefficient

# Just like the other load_problem() method except that the matrix is
# given in a standard column major ordered format (without gaps).
function load_problem(model::ClpModel,  num_cols::Integer, num_rows::Integer,
        start::Vector{CoinBigIndex}, index::Vector{Int32},
        value::Vector{Float64},
        col_lb::VecOrNothing, col_ub::VecOrNothing,
        obj::VecOrNothing,
        row_lb::VecOrNothing, row_ub::VecOrNothing)
    _jl__check_model(model)
    @clp_ccall loadProblem Void (Ptr{Void},Int32,Int32,Ptr{CoinBigIndex},Ptr{Int32},
    Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}) model.p num_cols num_rows start index value vec_or_null(Float64,col_lb,num_cols) vec_or_null(Float64,col_ub,num_cols) vec_or_null(Float64,obj,num_cols) vec_or_null(Float64,row_lb,num_rows) vec_or_null(Float64,row_ub,num_rows)
end

function load_problem(model::ClpModel,  constraint_matrix::AbstractMatrix, 
    col_lb::VecOrNothing, col_ub::VecOrNothing, 
    obj::VecOrNothing, row_lb::VecOrNothing,
    row_ub::VecOrNothing)
    mat = convert(SparseMatrixCSC{Float64,Int32},constraint_matrix)
    # We need to convert to zero-based, but
    # TODO: don't make extra copies of arrays
    # TODO: check dimensions match
    load_problem(model,mat.n, mat.m,mat.colptr-convert(Int32,1),
        mat.rowval-convert(Int32,1),mat.nzval,
        col_lb,col_ub,obj,row_lb,row_ub)
end

# Read quadratic part of the objective (the matrix part).
function load_quadratic_objective(model::ClpModel,
        num_cols::Integer, start::Vector{CoinBigIndex},
        col::Vector{Int32}, element::Vector{Float64})
    _jl__check_model(model)
    @clp_ccall loadQuadraticObjective Void (Ptr{Void},Int32,Ptr{CoinBigIndex},Ptr{Int32},Ptr{Float64}) model.p num_cols start column element
end

function load_quadratic_objective(model::ClpModel,
    hessian_matrix::SparseMatrixCSC{Float64,Int32})
    load_quadratic_objective(model, hessian_matrix.n, hessian_matrix.colptr-convert(Int32,1), hessian_matrix.rowval-convert(Int32,1),hessian_matrix.nzval)
end

# Read an mps file from the given filename.
function read_mps(model::ClpModel, mpsfile::AbstractString, keep_names::Bool, ignore_errors::Bool)
    _jl__check_model(model)
    _jl__check_file_is_readable(mpsfile)

    status = @clp_ccall readMps Int32 (Ptr{Void}, Ptr{UInt8}, Int32, Int32) model.p bytestring(mpsfile) keep_names ignore_errors
    if status != 0
        error("read_mps: error reading file $mpsfile")
    end
    return true
end
read_mps(model::ClpModel, mpsfile::AbstractString) = read_mps(model, mpsfile, true, false)

# Copy in integer information.
function copy_in_integer_information(model::ClpModel, information::Vector{UInt8})
    _jl__check_model(model)
    if length(information) != num_cols(model)
        error("Length of integer information must match number of columns")
    end
    @clp_ccall copyInIntegerInformation Void (Ptr{Void},Ptr{UInt8}) model.p information
end

# Drop integer information.
function delete_integer_information(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall deleteIntegerInformation Void (Ptr{Void},) model.p
end

# Resize rim part of model.
function resize(model::ClpModel, new_num_rows::Integer, new_num_cols::Integer)
    _jl__check_model(model)
    @clp_ccall resize Void (Ptr{Void},Int32,Int32) model.p new_num_rows new_num_cols
end

# Delete rows.
function delete_rows(model::ClpModel, which::Vector{Int32})
    _jl__check_model(model)
    @clp_ccall deleteRows Void (Ptr{Void},Int32,Ptr{Int32}) model.p length(which) which
end

# Add rows.
function add_rows(model::ClpModel, number::Integer, row_lower::Vector{Float64},
        row_upper::Vector{Float64},
        row_starts::Vector{Int32}, columns::Vector{Int32},
        elements::Vector{Float64})
    _jl__check_model(model)

    @clp_ccall addRows Void (Ptr{Void}, Int32, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}) model.p number row_lower row_upper row_starts columns elements
end

# Delete columns.
function delete_columns(model::ClpModel, which::Vector{Int32})
    _jl__check_model(model)
    @clp_ccall deleteColumns Void (Ptr{Void},Int32,Ptr{Int32}) model.p length(which) which
end

# Add columns.
function add_columns(model::ClpModel, number::Integer, column_lower::Vector{Float64},
        column_upper::Vector{Float64},
        objective::Vector{Float64},
        column_starts::Vector{Int32}, rows::Vector{Int32},
        elements::Vector{Float64})
    _jl__check_model(model)

    @clp_ccall addColumns Void (Ptr{Void}, Int32, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}) model.p number column_lower column_upper objective column_starts rows elements
    
end

function add_columns(model::ClpModel, column_lower::Vector{Float64},
        column_upper::Vector{Float64},
        objective::Vector{Float64},
        new_columns::SparseMatrixCSC{Float64,Int32})

    add_columns(model, new_columns.n, column_lower, column_upper, objective, new_columns.colptr-convert(Int32,1), new_columns.rowval-convert(Int32,1), new_columns.nzval)
end


# Change row lower bounds.
function chg_row_lower(model::ClpModel, row_lower::Vector{Float64})
    _jl__check_model(model)
    if length(row_lower) != get_num_rows(model)
        error("Array length must match number of rows in the model")
    end

    @clp_ccall chgRowLower Void (Ptr{Void}, Ptr{Float64}) model.p row_lower
end

# Change row upper bounds.
function chg_row_upper(model::ClpModel, row_upper::Vector{Float64})
    _jl__check_model(model)
    if length(row_upper) != get_num_rows(model)
        error("Array length must match number of rows in the model")
    end
    
    @clp_ccall chgRowUpper Void (Ptr{Void}, Ptr{Float64}) model.p row_upper
end

# Change column lower bounds.
function chg_column_lower(model::ClpModel, column_lower::Vector{Float64})
    _jl__check_model(model)
    if length(column_lower) != get_num_cols(model)
        error("Array length must match number of columns in the model")
    end

    @clp_ccall chgColumnLower Void (Ptr{Void}, Ptr{Float64}) model.p column_lower 
end

# Change column upper bounds.
function chg_column_upper(model::ClpModel, column_upper::Vector{Float64})
    _jl__check_model(model)
    if length(column_upper) != get_num_cols(model)
        error("Array length must match number of columns in the model")
    end

    @clp_ccall chgColumnUpper Void (Ptr{Void}, Ptr{Float64}) model.p column_upper 
end

# Change objective coefficients.
function chg_obj_coefficients(model::ClpModel, obj_in::Vector{Float64})
    _jl__check_model(model)
    if length(obj_in) != get_num_cols(model)
        error("Array length must match number of columns in the model")
    end
    
    @clp_ccall chgObjCoefficients Void (Ptr{Void},Ptr{Float64}) model.p obj_in
end

# Drops names - makes lengthnames 0 and names empty.
function drop_names(model::ClpModel)
    _jl__check_model(model)

    @clp_ccall dropNames Void (Ptr{Void},) model.p
end

# Copy in names.
# TODO: Need a nice way to accept a vector of Strings
function copy_names(model::ClpModel, row_names::Vector{Vector{UInt8}},
        columnNames::Vector{Vector{UInt8}})
    # TODO
    error("TODO")
    return
end

# Number of rows.
function get_num_rows(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall getNumRows Int32 (Ptr{Void},) model.p
end
number_rows = get_num_rows

# Number of columns.
function get_num_cols(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall getNumCols Int32 (Ptr{Void},) model.p
end
number_cols = get_num_cols

# Get primal tolerance.
function primal_tolerance(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall primalTolerance Float64 (Ptr{Void},) model.p
end

# Set primal tolerance to use.
function set_primal_tolerance(model::ClpModel, value::Float64)
    _jl__check_model(model)
    @clp_ccall setPrimalTolerance Void (Ptr{Void}, Float64) model.p value
    return
end

# Get dual tolerance.
function dual_tolerance(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall dualTolerance Float64 (Ptr{Void},) model.p
end

# Set dual tolerance to use.
function set_dual_tolerance(model::ClpModel, value::Float64)
    _jl__check_model(model)
    @clp_ccall setDualTolerance Void (Ptr{Void}, Float64) model.p value
    return
end

# Get dual objective limit.
function dual_objective_limit(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall dualObjectiveLimit Float64 (Ptr{Void},) model.p
end

# Set dual objective limit.
function set_dual_objective_limit(model::ClpModel, value::Float64)
    @clp_ccall setDualObjectiveLimit Void (Ptr{Void}, Float64) model.p value
    return
end

# Get objective offset.
function objective_offset(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall dualObjectiveLimit Float64 (Ptr{Void},) model.p
end

# Set objective offset.
function set_objective_offset(model::ClpModel, value::Float64)
    @clp_ccall setObjectiveOffset Void (Ptr{Void}, Float64) model.p value
    return
end

# Fill in array with problem name.
function problem_name(model::ClpModel)
    _jl__check_model(model)
    problem_name_p = pointer(Array(UInt8, 1000))
    @clp_ccall problemName Void (Ptr{Void}, Int32, Ptr{UInt8}) model.p 1000 problem_name_p
    return bytestring(problem_name_p)
end

# Set problem name.  Must have \0 at end.
function set_problem_name(model::ClpModel, name::String)
    _jl__check_model(model)
    @assert isascii(name)
    
    @clp_ccall setProblemName Void (Ptr{Void}, Int32, Ptr{UInt8}) model.p (length(name)+1) bytestring(name)
end

# Get number of iterations
function number_iterations(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall numberIterations Int32 (Ptr{Void},) model.p
end
get_iteration_count = number_iterations

# Set number of iterations
function set_number_iterations(model::ClpModel, iters::Integer)
    _jl__check_model(model)
    @clp_ccall setNumberIterations Void (Ptr{Void}, Int32) model.p iters
    return
end

# Get maximum number of iterations
function maximum_iterations(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall maximumIterations Int32 (Ptr{Void},) model.p
end

# Set maximum number of iterations
function set_maximum_iterations(model::ClpModel, max_iters::Integer)
    _jl__check_model(model)
    @clp_ccall setMaximumIterations Void (Ptr{Void}, Int32) model.p max_iters
    return
end

# Get maximum time in seconds (from when set is called)
function maximum_seconds(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall maximumSeconds Float64 (Ptr{Void},) model.p
end

# Set maximum time in seconds (from when set is called)
function set_maximum_seconds(model::ClpModel, max_secs::Real)
    _jl__check_model(model)
    @clp_ccall setMaximumSeconds Void (Ptr{Void}, Float64) model.p max_secs
    return
end

# Query whether maximum iterations or time were hit
function hit_maximum_iterations(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall hitMaximumIterations Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Get the status of the problem:
#   0 - optimal
#   1 - primal infeasible
#   2 - dual infeasible
#   3 - stopped on iterations etc
#   4 - stopped due to errors
function status(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall status Int32 (Ptr{Void},) model.p
end

# Set the status of the problem.
function set_problem_status(model::ClpModel, status::Integer)
    _jl__check_model(model)
    @clp_ccall setProblemStatus Void (Ptr{Void}, Int32) model.p status
    return
end

# Get the secondary status of problem - may get extended:
#   0 - none
#   1 - primal infeasible because dual limit reached
#   2 - scaled problem optimal - unscaled has primal infeasibilities
#   3 - scaled problem optimal - unscaled has dual infeasibilities
#   4 - scaled problem optimal - unscaled has both dual and primal infeasibilities
function secondary_status(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall secondaryStatus Int32 (Ptr{Void},) model.p
end

# Set the secondary status of the problem.
function set_secondary_status(model::ClpModel, status::Integer)
    _jl__check_model(model)
    @clp_ccall setSecondaryStatus Void (Ptr{Void}, Int32) model.p status
    return
end

# Get the direction of optimization:
#   1 - minimize
#   -1 - maximize
#   0 - ignore
# XXX: for some reason this returns a floating point (???)
function optimization_direction(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall optimizationDirection Float64 (Ptr{Void},) model.p
end

# Set the direction of optimization.
# XXX: for some reason this takes a floating point argument (???)
function set_optimization_direction(model::ClpModel, value::Real)
    _jl__check_model(model)
    @clp_ccall setOptimizationDirection Void (Ptr{Void}, Float64) model.p value
    return
end

# TODO: do we actually want to make a copy of the result?
# More efficient to not and let the sufficiently warned user do so if desired.
macro def_get_row_property(fname,clpname)
    quote
        function $(esc(fname))(model::ClpModel)
            _jl__check_model(model)
            row_property_p = @clp_ccall $clpname Ptr{Float64} (Ptr{Void},) model.p
            num_rows = convert(Int,get_num_rows(model))
            return copy(unsafe_wrap(Array,row_property_p,(num_rows,)))
        end
    end
end

macro def_get_col_property(fname,clpname)
    quote
        function $(esc(fname))(model::ClpModel)
            _jl__check_model(model)
            col_property_p = @clp_ccall $clpname Ptr{Float64} (Ptr{Void},) model.p
            num_cols = convert(Int,get_num_cols(model))
            return copy(unsafe_wrap(Array,col_property_p,(num_cols,)))
        end
    end
end

# Get primal row solution.
@def_get_row_property get_row_activity getRowActivity

primal_row_solution = get_row_activity

# Get primal column solution.
@def_get_col_property get_col_solution getColSolution
primal_column_solution = get_col_solution

function set_col_solution(model::ClpModel, input::Vector{Float64})
    _jl__check_model(model)

    @clp_ccall setColSolution Void (Ptr{Void},Ptr{Float64}) model.p input
end

# Get dual row solution.
@def_get_row_property get_row_price getRowPrice
dual_row_solution = get_row_price

# Get dual column solution (i.e. the reduced costs).
@def_get_col_property get_reduced_cost getReducedCost
dual_column_solution = get_reduced_cost

# Get row lower bounds.
@def_get_row_property get_row_lower getRowLower
row_lower = get_row_lower

# Get row upper bounds.
@def_get_row_property get_row_upper getRowUpper
row_upper = get_row_upper

@def_get_col_property get_obj_coefficients getObjCoefficients
objective = get_obj_coefficients

# Get column lower bounds.
@def_get_col_property get_col_lower getColLower
column_lower = get_col_lower

# Get column upper bounds.
@def_get_col_property get_col_upper getColUpper
column_upper = get_col_upper

# Get the number of elements in matrix.
function get_num_elements(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall getNumElements Int32 (Ptr{Void},) model.p
end

# Get column starts in matrix.
function get_vector_starts(model::ClpModel)
    _jl__check_model(model)
    vec_starts_p = @clp_ccall getVectorStarts Ptr{CoinBigIndex} (Ptr{Void},) model.p
    num_cols = Int(get_num_cols(model))
    return copy(unsafe_wrap(Array,vec_starts_p, (num_cols+1,)))
end

# Get row indices in matrix.
function get_indices(model::ClpModel)
    _jl__check_model(model)
    # getIndices returns an "int*", how do we know it's Int32??
    row_indices_p = @clp_ccall getIndices Ptr{Int32} (Ptr{Void},) model.p
    num_elts = Int(get_num_elements(model))
    return copy(unsafe_wrap(Array,row_indices_p,(num_elts,)))
end

function get_constraint_matrix(model::ClpModel)
    _jl__check_model(model)
    @assert CoinBigIndex == Int32
    # SparseMatrixCSC requires same integer type for colptr and rowval
    num_cols = convert(Int,get_num_cols(model))
    num_rows = convert(Int,get_num_rows(model))
    colptr = get_vector_starts(model) + convert(Int32,1)
    rowval = get_indices(model) + convert(CoinBigIndex,1)
    nzval = get_elements(model)

    return SparseMatrixCSC{Float64,CoinBigIndex}(num_rows,num_cols,colptr,rowval,nzval)

end

# Get column vector lengths in matrix.
@def_get_col_property get_vector_lengths getVectorLengths

# Get element values in matrix.
function get_elements(model::ClpModel)
    _jl__check_model(model)
    elements_p = @clp_ccall getElements Ptr{Float64} (Ptr{Void},) model.p
    num_elts = Int(get_num_elements(model))
    return copy(unsafe_wrap(Array,elements_p,(num_elts,)))
end

# Get objective value.
function get_obj_value(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall getObjValue Float64 (Ptr{Void},) model.p
end
objective_value = get_obj_value

# Get integer information.
function integer_information(model::ClpModel)
    # TODO
    error("TODO")
    return # vector of chars (?)
end

# Get an infeasibility ray (empty vector returned if none/wrong).
function infeasibility_ray(model::ClpModel)
    _jl__check_model(model)
    infeas_ray_p = @clp_ccall infeasibilityRay Ptr{Float64} (Ptr{Void},) model.p
    num_rows = convert(Int,get_num_rows(model))
    local infeas_ray::Vector{Float64}
    if infeas_ray_p != C_NULL
        infeas_ray = copy(unsafe_wrap(Array,infeas_ray_p,(num_rows,)))
        ccall(:free,Void,(Ptr{Void},),infeas_ray_p)
    else
        infeas_ray = Array(Float64,0)
    end
    return infeas_ray
end

# Get an unbounded ray (empty vector returned if none/wrong).
function unbounded_ray(model::ClpModel)
    _jl__check_model(model)
    unbd_ray_p = @clp_ccall unboundedRay Ptr{Float64} (Ptr{Void},) model.p
    num_cols = convert(Int,get_num_cols(model))
    local unbd_ray::Vector{Float64}
    if unbd_ray_p != C_NULL 
        unbd_ray = copy(unsafe_wrap(Array,unbd_ray_p,(num_cols,)))
        ccall(:free,Void,(Ptr{Void},),unbd_ray_p)
    else
        unbd_ray = Array(Float64,0)
    end
    return unbd_ray
end

# Query whether the status array exists (partly for OsiClp).
function statusExists(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall statusExists Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Return the address of the status array (char[num_rows+num_cols]).
# Status values are:
#   0 - free
#   1 - basic
#   2 - at upper
#   3 - at lower
#   4 - superbasic
#   (5 - fixed)
#   [???]
function status_array(model::ClpModel)
    # TODO
    error("TODO")
    return # pointer to unsigned chars (?)
end

# Copy in the status vector
# XXX copyin rather than copy_in, WTF?
# XXX input vector of unsigned chars (?)
function copyin_status(model::ClpModel, status_array::Vector{UInt8})
    # TODO
    error("TODO")
    return
end

# Get variable basis info.
# Note that we adjust from one-based indices
function get_column_status(model::ClpModel, sequence::Integer)
    _jl__check_model(model)
    _jl__check_col_is_valid(model, sequence)
    return @clp_ccall getColumnStatus Int32 (Ptr{Void},Int32) model.p (sequence-1)
end

# Get row basis info.
function get_row_status(model::ClpModel, sequence::Integer)
    _jl__check_model(model)
    _jl__check_row_is_valid(model, sequence)
    return @clp_ccall getRowStatus Int32 (Ptr{Void},Int32) model.p (sequence-1)
end

# Set variable basis info (and value if at bound).
function set_column_status(model::ClpModel, sequence::Integer, value::Integer)
    _jl__check_model(model)
    _jl__check_col_is_valid(model, sequence)
    @clp_ccall setColumnStatus Void (Ptr{Void},Int32,Int32) model.p (sequence-1) value
end

# Set row basis info (and value if at bound).
function set_row_status(model::ClpModel, sequence::Integer, value::Integer)
    _jl__check_model(model)
    _jl__check_row_is_valid(model, sequence)
    @clp_ccall setRowStatus Void (Ptr{Void},Int32,Int32) model.p (sequence-1) value
end

# No reason to expose these
# Set user pointer (used for generic purposes)
#function get_user_pointer(model::ClpModel)
# Get user pointer (used for generic purposes)
#function set_user_pointer(model::ClpModel, pointer::Ptr{Void})

# Pass in c-pointer to callback function. (from cfunction())
# See Coin_C_defines.h for signature
# TODO: test this
function register_call_back(model::ClpModel, callback::Ptr{Void})
    _jl__check_model(model)
    @clp_ccall registerCallBack Void (Ptr{Void},Ptr{Void}) model.p callback
end

# Unset callback function.
function clear_call_back(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall clearCallBack Void (Ptr{Void},) model.p
    return
end

# Get amount of print out:
#   0 - none
#   1 - just final
#   2 - just factorizations
#   3 - as 2 plus a bit more
#   4 - verbose
#   above that: 8,16,32 etc just for selective debug.
function log_level(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall logLevel Int32 (Ptr{Void},) model.p
end

# Set amount of print out.
function set_log_level(model::ClpModel, value::Integer)
    _jl__check_model(model)
    @clp_ccall setLogLevel Void (Ptr{Void}, Int32) model.p value
    return
end

# Get length of names (0 means no names).
function length_names(model::ClpModel)
    _jl__check_model(model)
    return @clp_ccall lengthNames Int32 (Ptr{Void},) model.p
end

# Return an array with a row name.
# XXX deviation from Clp API syntax, which fills a user-provided
# array of chars:
#   void Clp_rowName(Clp_Simplex * model, int iRow, char * name);
function row_name(model::ClpModel, row::Integer)
    _jl__check_model(model)
    _jl__check_row_is_valid(model, row)
    size = @clp_ccall lengthNames Int32 (Ptr{Void},) model.p
    row_name_p = pointer(Array(UInt8, size+1))
    @clp_ccall rowName Void (Ptr{Void}, Int32, Ptr{UInt8}) model.p (row-1) row_name_p
    return bytestring(row_name_p)
end

# Return an array with a column name.
# XXX deviation from Clp API syntax, which fills a user-provided
# array of chars:
#   void Clp_columnName(Clp_Simplex * model, int iColumn, char * name);
function column_name(model::ClpModel, col::Integer)
    _jl__check_model(model)
    _jl__check_col_is_valid(model, col)
    size = @clp_ccall lengthNames Int32 (Ptr{Void},) model.p
    col_name_p = pointer(Array(UInt8,size+1))
    @clp_ccall columnName Void (Ptr{Void}, Int32, Ptr{UInt8}) model.p (col-1) col_name_p
    return bytestring(col_name_p)
end

# General solve algorithm which can do presolve.
function initial_solve(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall initialSolve Int32 (Ptr{Void},) model.p
end

function initial_solve_with_options(model::ClpModel,solve::ClpSolve)
    _jl__check_model(model)
    _jl__check_solve(solve)
    @clp_ccall initialSolveWithOptions Int32 (Ptr{Void},Ptr{Void}) model.p solve.p
end

initial_solve(model::ClpModel, solve::ClpSolve) = initial_solve_with_options(model,solve)

# Dual initial solve.
function initial_dual_solve(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall initialDualSolve Int32 (Ptr{Void},) model.p
end

# Primal initial solve.
function initial_primal_solve(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall initialPrimalSolve Int32 (Ptr{Void},) model.p
end

# Barrier initial solve.
function initial_barrier_solve(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall initialBarrierSolve Int32 (Ptr{Void},) model.p
end

# Barrier initial solve, no crossover.
function initial_barrier_no_cross_solve(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall initialBarrierNoCrossSolve Int32 (Ptr{Void},) model.p
end

# Dual algorithm.
function dual(model::ClpModel, ifValuesPass::Integer)
    _jl__check_model(model)
    @clp_ccall dual Int32 (Ptr{Void},Int32) model.p ifValuesPass
end

# Primal algorithm.
function primal(model::ClpModel, ifValuesPass::Integer)
    _jl__check_model(model)
    @clp_ccall primal Int32 (Ptr{Void},Int32) model.p ifValuesPass
end

# Solve the problem with the idiot code.
function idiot(model::ClpModel, tryhard::Integer)
    _jl__check_model(model)
    @clp_ccall idiot Int32 (Ptr{Void},Int32) model.p tryhard
end

# Set or unset scaling:
#  0 - off
#  1 - equilibrium
#  2 - geometric
#  3 - auto
#  4 - dynamic (later)
function scaling(model::ClpModel, mode::Integer)
    _jl__check_model(model)
    if !(0 <= mode <= 4)
        error("Invalid clp scaling mode $mode (must be between 0 and 4)")
    end
    @clp_ccall scaling Void (Ptr{Void}, Int32) model.p mode
    return
end

# Get scaling flag.
function scaling_flag(model::ClpModel)
    _jl__check_model(model)
    return @clp_ccall scalingFlag Int32 (Ptr{Void},) model.p
end

# Crash (at present just aimed at dual); returns:
#   -2 if dual preferred and crash basis created
#   -1 if dual preferred and all slack basis preferred
#   0 if basis going in was not all slack
#   1 if primal preferred and all slack basis preferred
#   2 if primal preferred and crash basis created.
#
#   If gap between bounds <="gap" variables can be flipped
#
#   If "pivot" is:
#     0 - No pivoting (so will just be choice of algorithm)
#     1 - Simple pivoting e.g. gub
#     2 - Mini iterations
function crash(model::ClpModel, gap::Float64, pivot::Int32)
    _jl__check_model(model)
    if !(0 <= pivot <= 2)
        error("Invalid clp crash pivot value $pivot (must be between 0 and 2)")
    end
    return @clp_ccall crash Int32 (Ptr{Void}, Float64, Int32) model.p gap pivot
end

# Query whether the problem is primal feasible.
function primal_feasible(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall primalFeasible Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Query whether the problem is dual feasible.
function dual_feasible(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall dualFeasible Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Get dual bound.
function dual_bound(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall dualBound Float64 (Ptr{Void},) model.p
end

# Set dual bound.
function set_dual_bound(model::ClpModel, value::Real)
    _jl__check_model(model)
    @clp_ccall setDualBound Void (Ptr{Void}, Float64) model.p value
    return
end

# Get infeasibility cost.
function infeasibility_cost(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall infeasibilityCost Float64 (Ptr{Void},) model.p
end

# Set infeasibility cost.
function set_infeasibility_cost(model::ClpModel, value::Real)
    _jl__check_model(model)
    @clp_ccall setInfeasibilityCost Void (Ptr{Void}, Float64) model.p value
    return
end

# Get Perturbation. Values are:
#   50  - switch on perturbation
#   100 - auto perturb if takes too long (1.0e-6 largest nonzero)
#   101 - we are perturbed
#   102 - don't try perturbing again
# The default is 100; others are for playing.
function perturbation(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall perturbation Int32 (Ptr{Void},) model.p
end

# Set Perturbation.
function set_perturbation(model::ClpModel, value::Integer)
    _jl__check_model(model)
    if !(value == 50 || value == 100 || value == 101 || value == 102)
        error("Invalid clp perturbation value: $value (must be one of 50,100,101,102)")
    end
    @clp_ccall setPerturbation Void (Ptr{Void}, Int32) model.p value
    return
end

# Get current (or last) algorithm.
function algorithm(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall algorithm Int32 (Ptr{Void},) model.p
end

# Set algorithm.
function set_algorithm(model::ClpModel, value::Integer)
    _jl__check_model(model)
    # XXX which values of the algorithm are valid ???
    @clp_ccall setAlgorithm Void (Ptr{Void}, Int32) model.p value
    return
end

# Get the sum of dual infeasibilities.
function sum_dual_infeasibilities(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall sumDualInfeasibilities Float64 (Ptr{Void},) model.p
end

# Get the number of dual infeasibilities.
function number_dual_infeasibilities(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall numberDualInfeasibilities Int32 (Ptr{Void},) model.p
end

# Get the sum of primal infeasibilities.
function sum_primal_infeasibilities(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall sumPrimalInfeasibilities Float64 (Ptr{Void},) model.p
end

# Get the number of primal infeasibilities.
function number_primal_infeasibilities(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall numberPrimalInfeasibilities Int32 (Ptr{Void},) model.p
end

# Save model to file, returns 0 if success.  This is designed for
# use outside algorithms so does not save iterating arrays etc.
# It does not save any messaging information.
# Does not save scaling values.
# It does not know about all types of virtual functions.
function save_model(model::ClpModel, file_name::String)
    _jl__check_model(model)
    @assert isascii(file_name)
    
    @clp_ccall saveModel Int32 (Ptr{Void},Ptr{UInt8}) model.p bytestring(file_name)
end

# Restore model from file, returns 0 if success,
# deletes current model.
function restore_model(model::ClpModel, file_name::String)
    _jl__check_model(model)
    @assert isascii(file_name)

    @clp_ccall restoreModel Int32 (Ptr{Void},Ptr{UInt8}) model.p bytestring(file_name)
end

# Just check solution (for external use) - sets sum of
# infeasibilities etc.
function check_solution(model::ClpModel);
    _jl__check_model(model)
    @clp_ccall checkSolution Void (Ptr{Void},) model.p
    return
end

# Query whether there are numerical difficulties.
function is_abandoned(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall isAbandoned Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Query whether optimality is proven.
function is_proven_optimal(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall isProvenOptimal Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Query whether primal infeasiblity is proven.
function is_proven_primal_infeasible(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall isProvenPrimalInfeasible Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Query whether dual infeasiblity is proven.
function is_proven_dual_infeasible(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall isProvenDualInfeasible Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Query whether the given primal objective limit is reached.
function is_primal_objective_limit_reached(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall isPrimalObjectiveLimitReached Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Query whether the given dual objective limit is reached.
function is_dual_objective_limit_reached(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall isDualObjectiveLimitReached Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Query whether the iteration limit is reached
function is_iteration_limit_reached(model::ClpModel)
    _jl__check_model(model)
    ret = @clp_ccall isIterationLimitReached Int32 (Ptr{Void},) model.p
    return ret != 0
end

# Get the direction of optimization:
#   1  - minimize
#   -1 - maximize
#   0  - ignore
function get_obj_sense(model::ClpModel)
    _jl__check_model(model)
    @clp_ccall getObjSense Float64 (Ptr{Void},) model.p
end

# Set the direction of optimization.
function set_obj_sense(model::ClpModel, objsen::Real)
    _jl__check_model(model)
    if !(objsen == -1 || objsen == 0 || objsen == 1)
        error("Invalid clp objsense $objsense (should be -1,0 or 1)")
    end
    @clp_ccall setObjSense Void (Ptr{Void}, Float64) model.p objsen
    return
end

# Print model for debugging purposes
function print_model(model::ClpModel, prefix::String)
    _jl__check_model(model)
    @assert isascii(prefix)

    @clp_ccall printModel Void (Ptr{Void},Ptr{UInt8}) model.p bytestring(prefix)
end

function get_small_element_value(model::ClpModel)
    _jl__check_model(model)

    @clp_ccall getSmallElementValue Float64 (Ptr{Void},) model.p
end

function set_small_element_value(model::ClpModel,value::Float64)
    _jl__check_model(model)

    @clp_ccall setSmallElementValue Void (Ptr{Void},Float64) model.p value
end

# ClpSolve functions

function set_special_option(solve::ClpSolve, which::Integer, value::Integer, extraInfo::Integer)
    _jl__check_solve(solve)

    @clpsolve_ccall setSpecialOption Void (Ptr{Void},Int32,Int32,Int32) solve.p which value extraInfo
end

function get_special_option(solve::ClpSolve, which::Integer)
    _jl__check_solve(solve)

    @clpsolve_ccall getSpecialOption Int32 (Ptr{Void},Int32) solve.p which
end

macro def_get_int_property(fname,clpname)
    quote
        function $(esc(fname))(solve::ClpSolve)
            _jl__check_solve(solve)
            @clpsolve_ccall $clpname Int32 (Ptr{Void},) solve.p
        end
    end
end

macro def_set_int_property(fname,clpname)
    quote
        function $(esc(fname))(solve::ClpSolve, val::Integer)
            _jl__check_solve(solve)
            @clpsolve_ccall $clpname Void (Ptr{Void},Int32) solve.p val
        end
    end
end

# 0 - presolve on
# 1 - presolve off
# 2 - presolve number
# 3 - presolve number cost
function set_presolve_type(solve::ClpSolve, value::Integer)
    _jl__check_solve(solve)

    @clpsolve_ccall setPresolveType Void (Ptr{Void},Int32,Int32) solve.p value -1
end

@def_get_int_property get_presolve_type getPresolveType

# 0 - dual simplex
# 1 - primal simplex
# 2 - primal or sprint
# 3 - barrier
# 4 - barrier no crossover
# 5 - automatic
function set_solve_type(solve::ClpSolve, value::Integer)
    _jl__check_solve(solve)

    @clpsolve_ccall setSolveType Void (Ptr{Void},Int32,Int32) solve.p value -1
end

@def_get_int_property get_solve_type getSolveType

@def_get_int_property get_presolve_passes getPresolvePasses

function get_extra_info(solve::ClpSolve, which::Integer)
    _jl__check_solve(solve)

    @clpsolve_ccall getExtraInfo Int32 (Ptr{Void},Int32) solve.p which
end

@def_set_int_property set_infeasible_return setInfeasibleReturn
@def_get_int_property infeasible_return infeasibleReturn

@def_get_int_property do_dual doDual
@def_set_int_property set_do_dual setDoDual

@def_get_int_property do_singleton doSingleton
@def_set_int_property set_do_singleton setDoSingleton

@def_get_int_property do_doubleton doDoubleton
@def_set_int_property set_do_doubleton setDoDoubleton

@def_get_int_property do_tripleton doTripleton
@def_set_int_property set_do_tripleton setDoTripleton

@def_get_int_property do_tighten doTighten
@def_set_int_property set_do_tighten setDoTighten

@def_get_int_property do_forcing doForcing
@def_set_int_property set_do_forcing setDoForcing

@def_get_int_property do_implied_free doImpliedFree
@def_set_int_property set_do_implied_free setDoImpliedFree

@def_get_int_property do_dupcol doDupcol
@def_set_int_property set_do_dupcol setDoDupcol

@def_get_int_property do_duprow doDuprow
@def_set_int_property set_do_duprow setDoDuprow

@def_get_int_property do_singleton_column doSingletonColumn
@def_set_int_property set_do_singleton_column setDoSingletonColumn

@def_get_int_property presolve_actions presolveActions
@def_set_int_property set_presolve_actions setPresolveActions

@def_get_int_property substitution substitution
@def_set_int_property set_substitution setSubstitution

end
