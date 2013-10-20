
module ClpMathProgSolverInterface
using Clp.ClpCInterface

require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
importall MathProgSolverInterface

export ClpMathProgModel,
    ClpSolver,
    loadproblem,
    writeproblem,
    getvarLB,
    setvarLB,
    getvarLB,
    setvarLB,
    getconstrLB,
    setconstrLB,
    getconstrUB,
    setconstrUB,
    getobj,
    setobj,
    addvar,
    addconstr,
    updatemodel,
    setsense,
    getsense,
    numvar,
    numconstr,
    optimize,
    status,
    getobjval,
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver


type ClpMathProgModel <: AbstractMathProgModel
    inner::ClpModel
end

immutable ClpSolver <: AbstractMathProgSolver
    options 
end
ClpSolver(;kwargs...) = ClpSolver(kwargs)

function ClpMathProgModel(;kwargs...)
    if length(kwargs) != 0
        warn("ClpMathProgModelInterface does not yet support options")
    end
    m = ClpMathProgModel(ClpModel())
    set_log_level(m.inner,0)
    return m
end

model(s::ClpSolver) = ClpMathProgModel(;s.options...)


function loadproblem(m::ClpMathProgModel, filename::String)
    if ends_with(filename,".mps") || ends_with(filename,".mps.gz")
       read_mps(m.inner,filename)
    else
       error("unrecognized input format extension in $filename")
    end
end   


loadproblem(m::ClpMathProgModel, A, collb, colub, obj, rowlb, rowub) = 
    load_problem(m.inner,A,collb,colub,obj,rowlb,rowub)



#writeproblem(m, filename::String)

getvarLB(m::ClpMathProgModel) = get_col_lower(m.inner)
setvarLB(m::ClpMathProgModel, collb) = chg_column_lower(m.inner, collb)

getvarUB(m::ClpMathProgModel) = get_col_upper(m.inner)
setvarUB(m::ClpMathProgModel, colub) = chg_column_upper(m.inner, colub)

getconstrLB(m::ClpMathProgModel) = get_row_lower(m.inner)
setconstrLB(m::ClpMathProgModel, rowlb) = chg_row_lower(m.inner, rowlb)

getconstrUB(m::ClpMathProgModel) = get_row_upper(m.inner)
setconstrUB(m::ClpMathProgModel, rowub) = chg_row_upper(m.inner, rowub)

getobj(m::ClpMathProgModel) = get_obj_coefficients(m.inner)
setobj(m::ClpMathProgModel, obj) = chg_obj_coefficients(m.inner, obj)

function addvar(m::ClpMathProgModel, rowidx, rowcoef, collb, colub, objcoef)
    @assert length(rowidx) == length(rowcoef)
    colstarts = Int32[0, length(rowcoef)]
    rows = Int32[ i - 1 for i in rowidx ]
    add_columns(m.inner, 1, Float64[collb], Float64[colub], Float64[objcoef],
       colstarts, rows, convert(Vector{Float64},rowcoef))
end

function addconstr(m::ClpMathProgModel, colidx, colcoef, rowlb, rowub)
    @assert length(colidx) == length(colcoef)
    rowstarts = Int32[0, length(colcoef)]
    cols = Int32[ i - 1 for i in colidx ]
    add_rows(m.inner, 1, Float64[rowlb], Float64[rowub], rowstarts, cols, convert(Vector{Float64}, colcoef))
end

updatemodel(m::ClpMathProgModel) = nothing

function setsense(m::ClpMathProgModel,sense)
    if sense == :Min
        set_obj_sense(m.inner, 1.0)
    elseif sense == :Max
        set_obj_sense(m.inner, -1.0)
    else
        error("Unrecognized objective sense $sense")
    end
end

function getsense(m::ClpMathProgModel)
    s = get_obj_sense(m.inner)
    if s == 1.0
        return :Min
    elseif s == -1.0
        return :Max
    else
        error("Internal library error")
    end
end

numvar(m::ClpMathProgModel) = get_num_cols(m.inner) 
numconstr(m::ClpMathProgModel) = get_num_rows(m.inner)

optimize(m::ClpMathProgModel) = initial_solve(m.inner)

function status(m::ClpMathProgModel)
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

getobjval(m::ClpMathProgModel) = objective_value(m.inner)

getsolution(m::ClpMathProgModel) = primal_column_solution(m.inner) 
getconstrsolution(m::ClpMathProgModel) = primal_row_solution(m.inner)
getreducedcosts(m::ClpMathProgModel) = dual_column_solution(m.inner)

getconstrduals(m::ClpMathProgModel) = dual_row_solution(m.inner)


getrawsolver(m::ClpMathProgModel) = m.inner

end
