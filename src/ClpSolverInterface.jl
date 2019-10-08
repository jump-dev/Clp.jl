module ClpMathProgSolverInterface
using Clp.ClpCInterface
using LinearAlgebra
using SparseArrays

import MathProgBase
const MPB = MathProgBase

export ClpMathProgModel, ClpSolver

mutable struct ClpMathProgModel <: MPB.AbstractLinearQuadraticModel
    inner::ClpModel
    solveroptions::ClpSolve
end

struct ClpSolver <: MPB.AbstractMathProgSolver
    options
end
ClpSolver(;kwargs...) = ClpSolver(kwargs)

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

function setoption(m::ClpMathProgModel, name::Symbol, value)
    if haskey(optionmap, name)
        optionmap[name](m.inner,value)
    elseif haskey(solveoptionmap, name)
        solveoptionmap[name](m.solveroptions,value)
    else
        error("Unrecognized option: $name")
    end
end



function ClpMathProgModel(;kwargs...)
    m = ClpMathProgModel(ClpModel(),ClpSolve())
    set_log_level(m.inner,0) # disable output by default
    for (name, value) in kwargs
        setoption(m, name, value)
    end
    return m
end

MPB.LinearQuadraticModel(s::ClpSolver) = ClpMathProgModel(;s.options...)

MPB.ConicModel(s::ClpSolver) = MPB.LPQPtoConicBridge(MPB.LinearQuadraticModel(s))
MPB.supportedcones(s::ClpSolver) = [:Free,:Zero,:NonNeg,:NonPos]


function MPB.loadproblem!(m::ClpMathProgModel, filename::AbstractString)
    if endswith(filename,".mps") || endswith(filename,".mps.gz")
       read_mps(m.inner,filename)
    else
       error("unrecognized input format extension in $filename")
    end
end


function MPB.loadproblem!(m::ClpMathProgModel, A, collb, colub, obj, rowlb,
                          rowub, sense)
    load_problem(m.inner,A,collb,colub,obj,rowlb,rowub)
    MPB.setsense!(m, sense)
end



#writeproblem(m, filename::String)

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

MPB.getvarLB(m::ClpMathProgModel) = replaceInf(get_col_lower(m.inner))
MPB.setvarLB!(m::ClpMathProgModel, collb) = chg_column_lower(m.inner, collb)

MPB.getvarUB(m::ClpMathProgModel) = replaceInf(get_col_upper(m.inner))
MPB.setvarUB!(m::ClpMathProgModel, colub) = chg_column_upper(m.inner, colub)

MPB.getconstrLB(m::ClpMathProgModel) = replaceInf(get_row_lower(m.inner))
MPB.setconstrLB!(m::ClpMathProgModel, rowlb) = chg_row_lower(m.inner, rowlb)

MPB.getconstrUB(m::ClpMathProgModel) = replaceInf(get_row_upper(m.inner))
MPB.setconstrUB!(m::ClpMathProgModel, rowub) = chg_row_upper(m.inner, rowub)

MPB.getobj(m::ClpMathProgModel) = get_obj_coefficients(m.inner)
MPB.setobj!(m::ClpMathProgModel, obj) = chg_obj_coefficients(m.inner, obj)

function MPB.getconstrmatrix(m::ClpMathProgModel)
    A = get_constraint_matrix(m.inner)
    return convert(SparseMatrixCSC{Float64,Int},A)
end

function MPB.addvar!(m::ClpMathProgModel, rowidx, rowcoef, collb, colub, objcoef)
    @assert length(rowidx) == length(rowcoef)
    colstarts = Int32[0, length(rowcoef)]
    rows = Int32[ i - 1 for i in rowidx ]
    add_columns(m.inner, 1, Float64[collb], Float64[colub], Float64[objcoef],
       colstarts, rows, convert(Vector{Float64},rowcoef))
end

function MPB.addconstr!(m::ClpMathProgModel, colidx, colcoef, rowlb, rowub)
    @assert length(colidx) == length(colcoef)
    rowstarts = Int32[0, length(colcoef)]
    cols = Int32[ i - 1 for i in colidx ]
    add_rows(m.inner, 1, Float64[rowlb], Float64[rowub], rowstarts, cols, convert(Vector{Float64}, colcoef))
end

function MPB.setsense!(m::ClpMathProgModel,sense)
    if sense == :Min
        set_obj_sense(m.inner, 1.0)
    elseif sense == :Max
        set_obj_sense(m.inner, -1.0)
    else
        error("Unrecognized objective sense $sense")
    end
end

function MPB.getsense(m::ClpMathProgModel)
    s = get_obj_sense(m.inner)
    if s == 1.0
        return :Min
    elseif s == -1.0
        return :Max
    else
        error("Internal library error")
    end
end

MPB.numvar(m::ClpMathProgModel) = get_num_cols(m.inner)
MPB.numconstr(m::ClpMathProgModel) = get_num_rows(m.inner)

MPB.optimize!(m::ClpMathProgModel) = initial_solve_with_options(m.inner,m.solveroptions)

function MPB.status(m::ClpMathProgModel)
   s = ClpCInterface.status(m.inner)
   if s == 0
       return :Optimal
   elseif s == 1
       return :Infeasible
   elseif s == 2
       return :Unbounded
   elseif s == 3
       return :UserLimit
   elseif s == 4
       return :Error
   else
       error("Internal library error")
   end
end

MPB.getobjval(m::ClpMathProgModel) = objective_value(m.inner)

MPB.getsolution(m::ClpMathProgModel) = primal_column_solution(m.inner)
MPB.getconstrsolution(m::ClpMathProgModel) = primal_row_solution(m.inner)
MPB.getreducedcosts(m::ClpMathProgModel) = dual_column_solution(m.inner)

MPB.getconstrduals(m::ClpMathProgModel) = dual_row_solution(m.inner)

function MPB.getinfeasibilityray(m::ClpMathProgModel)
    return LinearAlgebra.rmul!(infeasibility_ray(m.inner),-1.0)
end
MPB.getunboundedray(m::ClpMathProgModel) = unbounded_ray(m.inner)

const statmap = Dict(zip([ 0x00,  0x01,            0x02,            0x03,       0x04,  0x05],
                     [:Free,:Basic,:NonbasicAtUpper,:NonbasicAtLower,:Superbasic,:Fixed]))
function MPB.getbasis(m::ClpMathProgModel)
    num_cols = MPB.numvar(m)
    num_rows = MPB.numconstr(m)
    cbasis = Array{Symbol}(undef, num_cols)
    rbasis = Array{Symbol}(undef, num_rows)
    for i in 1:num_cols
        val = get_column_status(m.inner, i)
        cbasis[i] = statmap[val]
    end
    for i in 1:num_rows
        val = get_row_status(m.inner, i)
        rbasis[i] = statmap[val]
    end
    return cbasis, rbasis
end

MPB.getvartype(m::ClpMathProgModel) = fill(:Cont, get_num_cols(m.inner))
function MPB.setvartype!(m::ClpMathProgModel, typ::Vector{Symbol})
    all(x->isequal(x,:Cont), typ) || error("Clp does not support integer variables")
    return
end

MPB.getrawsolver(m::ClpMathProgModel) = m.inner

end
