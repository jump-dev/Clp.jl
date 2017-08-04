module ClpMathProgSolverInterface
using Clp.ClpCInterface

importall MathProgBase.SolverInterface

export ClpMathProgModel,
    ClpSolver,
    loadproblem!,
    writeproblem,
    getvarLB,
    setvarLB!,
    getvarLB,
    setvarLB!,
    getconstrLB,
    setconstrLB!,
    getconstrUB,
    setconstrUB!,
    getobj,
    setobj!,
    getconstrmatrix,
    addvar!,
    addconstr!,
    setsense!,
    getsense,
    numvar,
    numconstr,
    optimize!,
    status,
    getobjval,
    getsolution,
    getconstrsolution,
    getreducedcosts,
    getconstrduals,
    getrawsolver


type ClpMathProgModel <: AbstractLinearQuadraticModel
    inner::ClpModel
    solveroptions::ClpSolve
end

immutable ClpSolver <: AbstractMathProgSolver
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

LinearQuadraticModel(s::ClpSolver) = ClpMathProgModel(;s.options...)

ConicModel(s::ClpSolver) = LPQPtoConicBridge(LinearQuadraticModel(s))
supportedcones(s::ClpSolver) = [:Free,:Zero,:NonNeg,:NonPos]


function loadproblem!(m::ClpMathProgModel, filename::AbstractString)
    if endswith(filename,".mps") || endswith(filename,".mps.gz")
       read_mps(m.inner,filename)
    else
       error("unrecognized input format extension in $filename")
    end
end   


function loadproblem!(m::ClpMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
    load_problem(m.inner,A,collb,colub,obj,rowlb,rowub)
    setsense!(m, sense)
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

getvarLB(m::ClpMathProgModel) = replaceInf(get_col_lower(m.inner))
setvarLB!(m::ClpMathProgModel, collb) = chg_column_lower(m.inner, collb)

getvarUB(m::ClpMathProgModel) = replaceInf(get_col_upper(m.inner))
setvarUB!(m::ClpMathProgModel, colub) = chg_column_upper(m.inner, colub)

getconstrLB(m::ClpMathProgModel) = replaceInf(get_row_lower(m.inner))
setconstrLB!(m::ClpMathProgModel, rowlb) = chg_row_lower(m.inner, rowlb)

getconstrUB(m::ClpMathProgModel) = replaceInf(get_row_upper(m.inner))
setconstrUB!(m::ClpMathProgModel, rowub) = chg_row_upper(m.inner, rowub)

getobj(m::ClpMathProgModel) = get_obj_coefficients(m.inner)
setobj!(m::ClpMathProgModel, obj) = chg_obj_coefficients(m.inner, obj)

function getconstrmatrix(m::ClpMathProgModel)
    A = get_constraint_matrix(m.inner)
    return convert(SparseMatrixCSC{Float64,Int},A)
end

function addvar!(m::ClpMathProgModel, rowidx, rowcoef, collb, colub, objcoef)
    @assert length(rowidx) == length(rowcoef)
    colstarts = Int32[0, length(rowcoef)]
    rows = Int32[ i - 1 for i in rowidx ]
    add_columns(m.inner, 1, Float64[collb], Float64[colub], Float64[objcoef],
       colstarts, rows, convert(Vector{Float64},rowcoef))
end

function addconstr!(m::ClpMathProgModel, colidx, colcoef, rowlb, rowub)
    @assert length(colidx) == length(colcoef)
    rowstarts = Int32[0, length(colcoef)]
    cols = Int32[ i - 1 for i in colidx ]
    add_rows(m.inner, 1, Float64[rowlb], Float64[rowub], rowstarts, cols, convert(Vector{Float64}, colcoef))
end

function setsense!(m::ClpMathProgModel,sense)
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

optimize!(m::ClpMathProgModel) = initial_solve_with_options(m.inner,m.solveroptions)

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

getinfeasibilityray(m::ClpMathProgModel) = scale!(infeasibility_ray(m.inner),-1.0)
getunboundedray(m::ClpMathProgModel) = unbounded_ray(m.inner)

const statmap = Dict(zip([ 0x00,  0x01,            0x02,            0x03,       0x04,  0x05],
                     [:Free,:Basic,:NonbasicAtUpper,:NonbasicAtLower,:Superbasic,:Fixed]))
function getbasis(m::ClpMathProgModel)
    num_cols = numvar(m)
    num_rows = numconstr(m)
    cbasis = Array{Symbol}(num_cols)
    rbasis = Array{Symbol}(num_rows)
    for i in 1:num_cols
        val = get_column_status(m.inner, i)
        cbasis[i] = statmap[val]
    end
    for i in 1:num_rows
        val = get_row_status(m.inner, i)
        rbasis[i] = statmap[val]
    end
    return cbasis,rbasis
end

getvartype(m::ClpMathProgModel) = fill(:Cont, get_num_cols(m.inner))
function setvartype!(m::ClpMathProgModel, typ::Vector{Symbol})
    all(x->isequal(x,:Cont), typ) || error("Clp does not support integer variables")
    return nothing
end

getrawsolver(m::ClpMathProgModel) = m.inner

end

module ClpMOISolverInterface
using Clp.ClpCInterface

# importall MathProgBase.SolverInterface

using MathOptInterface
const MOI = MathOptInterface

export ClpMOISolver

# Solver

immutable ClpMOISolver <: MOI.AbstractSolver
    options
end
ClpMOISolver(;kwargs...) = ClpMOISolver(kwargs)

function MOI.supportsproblem(s::ClpMOISolver, objective_type, constraint_types::Vector)
    if objective_type != MOI.ScalarAffineFunction{Float64}
        return false
    end
    supported_constraint_types = [(MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}),
            (MOI.ScalarAffineFunction{Float64}, MOI.GreaterThan{Float64}),
            (MOI.ScalarAffineFunction{Float64}, MOI.Interval{Float64}),
            (MOI.SingleVariable, MOI.LessThan{Float64}),
            (MOI.SingleVariable, MOI.GreaterThan{Float64})]
    for v in constraint_types
        if !(v in supported_constraint_types)
            return false
        end
    end
    return true
end

MOI.getattribute(s::ClpMOISolver, ::MOI.SupportsDuals) = true

# Solver instance

const VR = MOI.VariableReference
const CR = MOI.ConstraintReference

type ClpMOIModel <: MOI.AbstractSolverInstance
    inner::ClpModel
    solveroptions::ClpSolve
    nvariables::UInt64
    vartolbref::Vector{CR}
    vartoubref::Vector{CR}
    ubreftovar::Dict{CR, Int32}
    lbreftovar::Dict{CR, Int32}
    nupperbounds::UInt64
    nlowerbounds::UInt64
    nconstraints::UInt64 # does take into account bounds
    # constraints::Vector{CR}
    nrows_greaterthan::UInt64
    nrows_lessthan::UInt64
    nrows_interval::UInt64
    nrows::UInt64
    rowtorangeref::Vector{CR}
    rangereftorow::Dict{CR, Int32}
    solution::Vector{Float64}
    rowsolution::Vector{Float64}
    dualcolsolution::Vector{Float64}
    dualrowsolution::Vector{Float64}
end
function ClpMOIModel(;kwargs...)
    m = ClpMOIModel(ClpModel(),ClpSolve(), 0, Vector{CR}(), Vector{CR}(), Dict{CR, Int32}(), Dict{CR, Int32}(),
            0, 0, 0, 0, 0, 0, 0, Vector{CR}(), Dict{CR, Int32}(),
            Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}() )
    set_log_level(m.inner,0) # disable output by default
    for (name, value) in kwargs
        setoption(m, name, value)
    end
    return m
end

col(m::ClpMOIModel,ref::VR) = ref.value - 1 # column index in solver
row(m::ClpMOIModel,ref::CR) = m.rangereftorow[ref] - 1 # row index in solver

MOI.SolverInstance(s::ClpMOISolver) = ClpMOIModel(;s.options...)

# Variables
MOI.getattribute(m::ClpMOIModel, ::MOI.NumberOfVariables) = m.nvariables

function MOI.addvariable!(m::ClpMOIModel)
    numberInColumn = 0
    rows = Vector{Int32}()
    elements = Vector{Float64}()
    add_column(m.inner, numberInColumn,  rows, elements, -Inf, +Inf, 0.0)
    m.nvariables += 1
    m.nconstraints += 2
    nv = m.nvariables
    nc = m.nconstraints
    vref = VR(nv)
    lbref = CR{MOI.SingleVariable, MOI.GreaterThan{Float64}}(nc - 1)
    ubref = CR{MOI.SingleVariable, MOI.LessThan{Float64}}(nc)
    push!(m.vartolbref, lbref)
    push!(m.vartoubref, ubref)
    m.lbreftovar[lbref] = nv
    m.ubreftovar[ubref] = nv
    return vref
end

function MOI.addvariables!(m::ClpMOIModel, n::Integer)
    return [MOI.addvariable!(m) for i in 1:n]
end

function MOI.delete!(m::ClpMOIModel, vr::VR)
    which = Vector{Int32}()
    push!(which, col(m,vr))
    delete_columns(m.inner, which)
end


# Objective
function MOI.setobjective!(m::ClpMOIModel, sense::MOI.OptimizationSense, f::MOI.ScalarAffineFunction{Float64})
    obj_in = zeros(Float64, m.nvariables)
    for i in 1:length(f.variables)
        index = col(m,f.variables[i]) + 1
        coeff = f.coefficients[i]
        obj_in[index] = coeff
    end
    chg_obj_coefficients(m.inner, obj_in)
    set_objective_offset(m.inner, f.constant)
    if sense == MOI.MinSense
        set_obj_sense(m.inner, 1.0)
    elseif sense == MOI.MaxSense
        set_obj_sense(m.inner, -1.0)
    else
        set_obj_sense(m.inner, 0.0)
    end
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.ObjectiveSense)
    s = get_obj_sense(m.inner)
    if s == 1.0
        return MOI.MinSense
    elseif s == -1.0
        return MOI.MaxSense
    else
        error("library get_obj_sense returned an incorrect value")
    end
end

# MOI.getattribute(m::ClpMOIModel, ::MOI.ObjectiveFunction)
# MOI.modifyobjective!(m::ClpMOIModel, change::MOI.AbstractFunctionModification)

# Constraint
function _addrangeref!(F,S,m::ClpMOIModel, haslessthan::Bool, hasgreaterthan::Bool)
    @assert haslessthan || hasgreaterthan
    m.nconstraints += 1
    m.nrows += 1
    if !haslessthan
        m.nrows_lessthan += 1
    elseif !hasgreaterthan
        m.nrows_greaterthan += 1
    else
        m.nrows_interval += 1
    end
    cref = CR{F,S}(m.nconstraints)
    push!(m.rowtorangeref, cref)
    m.rangereftorow[cref] = m.nrows
    return cref
end

function _addrangeinsolver!(m::ClpMOIModel, f::MOI.ScalarAffineFunction{Float64}, lower::Float64, upper::Float64)
    elements = Vector{Float64}()
    columns = Vector{Int32}()
    number_in_row = length(f.variables)
    for i in 1:number_in_row
        push!(columns, col(m,f.variables[i]))
        push!(elements, f.coefficients[i])
    end
    add_row(m.inner, number_in_row, columns, elements, lower - f.constant, upper - f.constant)
end

function _addrange!(m::ClpMOIModel, f::MOI.ScalarAffineFunction{Float64}, lower::Float64, upper::Float64,
        haslower::Bool, hasupper::Bool)
    _addrangeinsolver!(m, f, lower, upper)
    return _addrangeref!(MOI.ScalarAffineFunction{Float64},MOI.Interval{Float64},m, haslower, hasupper)
end

function MOI.addconstraint!(m::ClpMOIModel, f::MOI.ScalarAffineFunction{Float64}, s::MOI.Interval{Float64})
    return _addrange!(m, f, s.lower, s.upper, true, true)
end

function MOI.addconstraint!(m::ClpMOIModel, f::MOI.ScalarAffineFunction{Float64}, s::MOI.GreaterThan{Float64})
    return _addrange!(m, f, s.lower, Inf, true, false)
end

function MOI.addconstraint!(m::ClpMOIModel, f::MOI.ScalarAffineFunction{Float64}, s::MOI.LessThan{Float64})
    return _addrange!(m, f, -Inf, s.upper, false, true)
end

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

function MOI.addconstraint!(m::ClpMOIModel, f::MOI.SingleVariable, s::MOI.GreaterThan{<:Real})
    lowerbounds = replaceInf(get_col_lower(m.inner))
    lowerbounds[col(m, f.variable) + 1] = s.lower
    m.nlowerbounds += 1
    chg_column_lower(m.inner, lowerbounds)
    return m.vartolbref[col(m, f.variable) + 1]
end

function MOI.addconstraint!(m::ClpMOIModel, f::MOI.SingleVariable, s::MOI.LessThan{<:Real})
    upperbounds = replaceInf(get_col_upper(m.inner))
    upperbounds[col(m, f.variable) + 1] = s.upper
    m.nupperbounds += 1
    chg_column_upper(m.inner, upperbounds)
    return m.vartoubref[col(m, f.variable) + 1]
end

function MOI.getattribute(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.SingleVariable, MOI.LessThan{Float64}})
    return m.nupperbounds
end

function MOI.getattribute(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.SingleVariable, MOI.GreaterThan{Float64}})
    return m.nlowerbounds
end

function MOI.getattribute(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}})
    return m.nrows_lessthan
end

function MOI.getattribute(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.GreaterThan{Float64}})
    return m.nrows_greaterthan
end

function MOI.getattribute(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.Interval{Float64}})
    return m.nrows_interval
end

# MOI.delete!(m::ClpMOIModel, cr::CR)
# MOI.modifyconstraint!(m::ClpMOIModel, cr::CR, change)

# Optimization and additional attributes
function MOI.optimize!(m::ClpMOIModel)
    initial_solve_with_options(m.inner,m.solveroptions)
    m.solution = primal_column_solution(m.inner)
    m.rowsolution = primal_row_solution(m.inner)
    m.dualcolsolution = dual_column_solution(m.inner) .* get_obj_sense(m.inner)
    m.dualrowsolution = dual_row_solution(m.inner) .* get_obj_sense(m.inner)
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.TerminationStatus)
    s = ClpCInterface.status(m.inner)
    if s == 0
        return MOI.Success
    elseif s == 1
        return MOI.InfeasibleNoResult
    elseif s == 2
        return MOI.UnboundedNoResult
    elseif s == 3
        return MOI.OtherLimit
    elseif s == 4
        return MOI.OtherError
    else
        error("Internal library error")
    end
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.PrimalStatus)
    if primal_feasible(m.inner)
        return MOI.FeasiblePoint
    else
        return MOI.UnknownResultStatus
    end
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.DualStatus)
    if dual_feasible(m.inner)
        return MOI.FeasiblePoint
    else
        return MOI.UnknownResultStatus
    end
end

MOI.getattribute(m::ClpMOIModel, ::MOI.ObjectiveValue) = objective_value(m.inner)

function MOI.getattribute(m::ClpMOIModel, ::MOI.VariablePrimal, vref::VR)
    return m.solution[col(m,vref)+1]
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.VariablePrimal, vec::Vector{VR})
    return [m.solution[col(m,vref)+1] for vref in vec]
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintPrimal, cref::CR)
    return m.rowsolution[row(m,cref)+1]
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintPrimal, vec::Vector{CR})
    return [m.rowsolution[row(m,cref)+1] for cref in vec]
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintDual, cref::CR{MOI.ScalarAffineFunction{Float64},S}
        where S <: Union{MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64}})
    return m.dualrowsolution[row(m,cref)+1]
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintDual,
        cref::CR{MOI.SingleVariable, MOI.GreaterThan{Float64}})
    var = m.lbreftovar[cref]
    if m.solution[var] == get_col_lower(m.inner)[var] #this wouold be more efficient to get from instance
        return m.dualcolsolution[var]
    else
        return 0
    end
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintDual,
        cref::CR{MOI.SingleVariable, MOI.LessThan{Float64}})
    var = m.ubreftovar[cref]
    if m.solution[var] == get_col_upper(m.inner)[var] #this wouold be more efficient to get from instance
        return m.dualcolsolution[var]
    else
        return 0
    end
end

function MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintDual, vec::Vector{CR})
    return [MOI.getattribute(m, MOI.ConstraintDual(), cref) for cref in vec]
end

MOI.cangetattribute(m::ClpMOIModel, ::Union{MOI.TerminationStatus,
                                    MOI.PrimalStatus,
                                    MOI.DualStatus,
                                    MOI.NumberOfVariables,
                                    MOI.NumberOfConstraints,
                                    MOI.ObjectiveSense,
                                    MOI.ObjectiveValue
                                    }) = true

MOI.cangetattribute(m::ClpMOIModel, ::Union{MOI.ConstraintPrimal, MOI.VariablePrimal, MOI.ConstraintDual},
        ::Union{VR, CR}) = true

MOI.cangetattribute(m::ClpMOIModel, ::Union{MOI.ConstraintPrimal, MOI.VariablePrimal},
        ::Union{Vector{VR}, Vector{CR}}) = true

MOI.cangetattribute(m::ClpMOIModel, ::MOI.ListOfConstraints) = false
MOI.cangetattribute(m::ClpMOIModel, ::Union{MOI.ConstraintFunction,
                                                 MOI.ConstraintSet}, ref::Union{VR, CR}) = false
MOI.cangetattribute(m::ClpMOIModel, ::MOI.ObjectiveFunction) = false
# MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintFunction, cr::CR)
# MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintSet, cr::CR)

end
