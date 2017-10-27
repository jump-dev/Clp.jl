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

MOI.get(s::ClpMOISolver, ::MOI.SupportsDuals) = true

# Solver instance

const VR = MOI.VariableReference
const CR = MOI.ConstraintReference
const IGL = Union{MOI.Interval{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}} 

type ClpMOIModel <: MOI.AbstractSolverInstance
    inner::ClpModel
    solveroptions::ClpSolve
    nvariables::UInt64
    coltovref::Vector{VR} # col means the column idx in solver + 1
    vreftocol::Dict{VR, Int32}
    vartolbref::Vector{CR}
    vartoubref::Vector{CR}
    vartoassignref::Vector{CR}
    ubreftovar::Dict{CR, Int32}
    lbreftovar::Dict{CR, Int32}
    assignreftovar::Dict{CR, Int32}
    nupperbounds::UInt64
    nlowerbounds::UInt64
    nassigns::UInt64
    nconstraints::UInt64 # Counts all cref created used or not, in solver or not.
    # constraints::Vector{CR}
    nrows_greaterthan::UInt64
    nrows_lessthan::UInt64
    nrows_equalto::UInt64
    nrows_interval::UInt64
    nrows::UInt64
    rowtorangeref::Vector{CR} # row means the row idx in solver + 1
    rangereftorow::Dict{CR, Int32}
    rangereftofunc::Dict{CR, MOI.ScalarAffineFunction{Float64}}
    rangereftolb::Dict{CR, Float64}
    rangereftoub::Dict{CR, Float64}    
    solution::Vector{Float64}
    rowsolution::Vector{Float64}
    dualcolsolution::Vector{Float64}
    dualrowsolution::Vector{Float64}
    objfunc::Nullable{MOI.ScalarAffineFunction{Float64}}
end

function ClpMOIModel(;kwargs...)
    m = ClpMOIModel(ClpModel(),ClpSolve(), 0, Vector{VR}(), Dict{VR, Int32}(),
            Vector{CR}(), Vector{CR}(), Vector{CR}(), 
            Dict{CR, Int32}(), Dict{CR, Int32}(), Dict{CR, Int32}(),
            0, 0, 0, 0, 0, 0, 0, 0, 0, Vector{CR}(), Dict{CR, Int32}(),
            Dict{CR, MOI.ScalarAffineFunction{Float64}}(), 
            Dict{CR, Float64}(), Dict{CR, Float64}(),
            Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(),
            Nullable{MOI.ScalarAffineFunction{Float64}}())
    set_log_level(m.inner,0) # disable output by default
    for (name, value) in kwargs
        setoption(m, name, value)
    end
    return m
end

# TODO : rename the following to scol and srow 
# (indices in solver are shifted by one compared to those in ClpMOIModel)
col(m::ClpMOIModel,ref::VR) = m.vreftocol[ref] - 1 # column index in solver 
row(m::ClpMOIModel,ref::CR) = m.rangereftorow[ref] - 1 # row index in solver

MOI.SolverInstance(s::ClpMOISolver) = ClpMOIModel(;s.options...)

# Variables
MOI.get(m::ClpMOIModel, ::MOI.NumberOfVariables) = m.nvariables

function MOI.addvariable!(m::ClpMOIModel)
    numberInColumn = 0
    rows = Vector{Int32}()
    elements = Vector{Float64}()
    add_column(m.inner, numberInColumn,  rows, elements, -Inf, +Inf, 0.0)
    m.nvariables += 1    
    m.nconstraints += 3
    nv = m.nvariables
    nc = m.nconstraints
    vref = VR(nv)
    push!(m.coltovref, vref)
    m.vreftocol[vref] = nv
    lbref = CR{MOI.SingleVariable, MOI.GreaterThan{Float64}}(nc - 2)
    ubref = CR{MOI.SingleVariable, MOI.LessThan{Float64}}(nc-1)
    assignref = CR{MOI.SingleVariable, MOI.EqualTo{Float64}}(nc)
    push!(m.vartolbref, lbref)
    push!(m.vartoubref, ubref)
    push!(m.vartoassignref, assignref)
    m.lbreftovar[lbref] = nv
    m.ubreftovar[ubref] = nv
    m.assignreftovar[assignref] = nv
    return vref
end

function MOI.addvariables!(m::ClpMOIModel, n::Integer)
    return [MOI.addvariable!(m) for i in 1:n]
end

function MOI.delete!(m::ClpMOIModel, vr::VR)
    column = col(m,vr)      
    m.nvariables -= 1
    
    deleteat!(m.coltovref, column+1)    
    for k in keys(m.vreftocol)
        if m.vreftocol[k] > column+1
            m.vreftocol[k] -= 1
        end
    end
    
    delete!(m.vreftocol, vr)    
    delete!(m.ubreftovar, m.vartoubref[vr.value])
    delete!(m.lbreftovar, m.vartolbref[vr.value])
    delete!(m.assignreftovar, m.vartoassignref[vr.value])
    
    lb = replaceInf(get_col_lower(m.inner))[column + 1]
    ub = replaceInf(get_col_upper(m.inner))[column + 1]    
    if lb == ub
        m.nassigns -= 1
    else
        if lb != -Inf
            m.nlowerbounds -= 1
        end
        if ub != Inf
            m.nupperbounds -= 1
        end
    end
    
    delete_columns(m.inner, Vector{Int32}([column]))
end

function MOI.set!(m::ClpMOIModel, ::MOI.ObjectiveFunction, f::MOI.ScalarAffineFunction{Float64})
    obj_in = zeros(Float64, m.nvariables)
    for i in 1:length(f.variables)
        index = col(m,f.variables[i]) + 1
        coeff = f.coefficients[i]
        obj_in[index] = coeff
    end
    m.objfunc = Nullable{MOI.ScalarAffineFunction{Float64}}(
            MOI.ScalarAffineFunction{Float64}(copy(f.variables),copy(f.coefficients), f.constant) )
    chg_obj_coefficients(m.inner, obj_in)
    set_objective_offset(m.inner, f.constant)    
end

function MOI.set!(m::ClpMOIModel, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    if sense == MOI.MinSense
        set_obj_sense(m.inner, 1.0)
    elseif sense == MOI.MaxSense
        set_obj_sense(m.inner, -1.0)
    else
        set_obj_sense(m.inner, 0.0)
    end
end

function MOI.get(m::ClpMOIModel, ::MOI.ObjectiveSense)
    s = get_obj_sense(m.inner)
    if s == 1.0
        return MOI.MinSense
    elseif s == -1.0
        return MOI.MaxSense
    else
        error("library get_obj_sense returned an incorrect value")
    end
end

function MOI.modifyobjective!(m::ClpMOIModel, change::MOI.ScalarCoefficientChange{Float64})
    @assert !isnull(m.objfunc)
    f = get(m.objfunc)        
    
    # modify function              
    varidx = findfirst(f.variables, change.variable)
    if varidx != 0
        f.coefficients[varidx] = change.new_coefficient
    else        
        push!(f.variables, change.variable)
        push!(f.coefficients, change.new_coefficient)
    end        
    
    # update instance in solver
    MOI.set!(m, MOI.ObjectiveFunction(), f)
end

function MOI.modifyobjective!(m::ClpMOIModel, change::MOI.ScalarConstantChange{Float64})    
    @assert !isnull(m.objfunc)
    f = get(m.objfunc)        
    
    # modify function              
    f.constant = change.new_constant    
    
    # update instance in solver
    MOI.set!(m, MOI.ObjectiveFunction(), f)
end

# Constraint
function _addrangeref!(f::MOI.ScalarAffineFunction{Float64},m::ClpMOIModel,
        lower::Float64, upper::Float64, haslower::Bool, hasupper::Bool)    
    @assert hasupper || haslower
    S = MOI.Interval{Float64}
    m.nconstraints += 1
    m.nrows += 1
    if !hasupper
        m.nrows_greaterthan += 1
        S = MOI.GreaterThan{Float64}
    elseif !haslower
        m.nrows_lessthan += 1
        S = MOI.LessThan{Float64}
    elseif lower == upper
        S = MOI.EqualTo{Float64}
        m.nrows_equalto += 1
    else    
        m.nrows_interval += 1
    end
    cref = CR{MOI.ScalarAffineFunction{Float64},S}(m.nconstraints)
    push!(m.rowtorangeref, cref)
    m.rangereftorow[cref] = m.nrows
    m.rangereftolb[cref] = lower
    m.rangereftoub[cref] = upper
    m.rangereftofunc[cref] = MOI.ScalarAffineFunction{Float64}(copy(f.variables), copy(f.coefficients), f.constant)
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
    return _addrangeref!(f,m, lower, upper, haslower, hasupper)
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

function MOI.addconstraint!(m::ClpMOIModel, f::MOI.ScalarAffineFunction{Float64}, s::MOI.EqualTo{Float64})
    return _addrange!(m, f, s.value, s.value, true, true)
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

function MOI.addconstraint!(m::ClpMOIModel, f::MOI.SingleVariable, s::MOI.GreaterThan{Float64})
    lowerbounds = replaceInf(get_col_lower(m.inner))
    lowerbounds[col(m, f.variable) + 1] = s.lower    
    chg_column_lower(m.inner, lowerbounds)
    m.nlowerbounds += 1
    return m.vartolbref[col(m, f.variable) + 1]
end

function MOI.addconstraint!(m::ClpMOIModel, f::MOI.SingleVariable, s::MOI.LessThan{Float64})
    upperbounds = replaceInf(get_col_upper(m.inner))
    upperbounds[col(m, f.variable) + 1] = s.upper    
    chg_column_upper(m.inner, upperbounds)
    m.nupperbounds += 1
    return m.vartoubref[col(m, f.variable) + 1]
end

function MOI.addconstraint!(m::ClpMOIModel, f::MOI.SingleVariable, s::MOI.EqualTo{Float64})
    upperbounds = replaceInf(get_col_upper(m.inner))
    upperbounds[col(m, f.variable) + 1] = s.value
    chg_column_upper(m.inner, upperbounds)
    lowerbounds = replaceInf(get_col_lower(m.inner))
    lowerbounds[col(m, f.variable) + 1] = s.value
    chg_column_lower(m.inner, lowerbounds)
    m.nassigns += 1    
    return m.vartoassignref[col(m, f.variable) + 1]
end

function MOI.modifyconstraint!(m::ClpMOIModel, 
        cref::CR{MOI.SingleVariable, MOI.GreaterThan{Float64}}, s::MOI.GreaterThan{Float64})
        
    lowerbounds = replaceInf(get_col_lower(m.inner))
    lowerbounds[m.lbreftovar[cref]] = s.lower    
    chg_column_lower(m.inner, lowerbounds)    
end

function MOI.modifyconstraint!(m::ClpMOIModel, 
        cref::CR{MOI.SingleVariable, MOI.LessThan{Float64}}, s::MOI.LessThan{Float64})
        
    upperbounds = replaceInf(get_col_upper(m.inner))
    upperbounds[m.ubreftovar[cref]] = s.upper    
    chg_column_upper(m.inner, upperbounds)    
end

function _modifyconstraintinsolver(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64}, S}  where S <: IGL, 
        f::MOI.ScalarAffineFunction{Float64}, lower::Float64, upper::Float64)        
    
    # update instance in solver
    r = row(m,cref)    
    delete_rows(m.inner, Vector{Int32}([r]))               
    _addrangeinsolver!(m, f::MOI.ScalarAffineFunction{Float64}, lower, upper)
    
    # update row <-> cref mapping.
    deleteat!(m.rowtorangeref, findfirst(m.rowtorangeref, cref))
    push!(m.rowtorangeref, cref)     
    m.rangereftorow[cref] = m.nrows            
end
        

function MOI.modifyconstraint!(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64}, S}  where S <: IGL, 
        change::MOI.ScalarCoefficientChange{Float64})
        
    lower = m.rangereftolb[cref]
    upper = m.rangereftoub[cref]
    f = m.rangereftofunc[cref]
            
    # modify function          
    varidx = findfirst(f.variables, change.variable)
    if varidx != 0
        f.coefficients[varidx] = change.new_coefficient
    else        
        push!(f.variables, change.variable)
        push!(f.coefficients, change.new_coefficient)
    end 
    
    _modifyconstraintinsolver(m, cref, f, lower, upper)               
end 

function _updatedetailednrows(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64},MOI.Interval{Float64}})
    m.nrows_interval -= 1    
end

function _updatedetailednrows(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64},MOI.GreaterThan{Float64}})
    m.nrows_greaterthan -= 1    
end

function _updatedetailednrows(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}})    
    m.nrows_lessthan -= 1
end

function _modifyrangeboundsinsolver(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64}, S} where S <: IGL)
        
    lower = replaceInf(get_row_lower(m.inner))
    upper = replaceInf(get_row_upper(m.inner))
        
    lower[row(m, cref)+1] = m.rangereftolb[cref]
    upper[row(m, cref)+1] = m.rangereftoub[cref]
    
    chg_row_lower(m.inner, lower)
    chg_row_upper(m.inner, upper)
end

function MOI.modifyconstraint!(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64}, S} where S <: IGL, 
        s::MOI.GreaterThan{Float64})
        
    # modify bounds in interface  
    _updatedetailednrows(m, cref)    
    m.nrows_greaterthan += 1
    m.rangereftolb[cref] = s.lower
    m.rangereftoub[cref] = Inf 
    
    _modifyrangeboundsinsolver(m, cref)
end

function MOI.modifyconstraint!(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64}, S} where S <: IGL, 
        s::MOI.LessThan{Float64})
        
    # modify bounds in interface      
    _updatedetailednrows(m, cref)
    
    m.nrows_lessthan += 1
    m.rangereftolb[cref] = -Inf
    m.rangereftoub[cref] = s.upper
    
    _modifyrangeboundsinsolver(m, cref)
end

function MOI.modifyconstraint!(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64}, S} where S <: IGL, 
        s::MOI.Interval{Float64})
        
    # modify bounds in interface  
    _updatedetailednrows(m, cref)    
    m.nrows_interval += 1
    m.rangereftolb[cref] = s.lower
    m.rangereftoub[cref] = s.upper
    
    _modifyrangeboundsinsolver(m, cref)
end

function MOI.delete!(m::ClpMOIModel, 
        cref::CR{MOI.ScalarAffineFunction{Float64}, S} where S <: IGL)
    
    _updatedetailednrows(m,cref)
    r = row(m,cref)
    m.nrows -= 1                   
         
    deleteat!(m.rowtorangeref, r+1)    
    for k in keys(m.rangereftorow)
        if m.rangereftorow[k] > r+1
            m.rangereftorow[k] -= 1
        end
    end
    
    delete!(m.rangereftorow, cref)
    delete!(m.rangereftolb, cref)
    delete!(m.rangereftoub, cref)
    delete!(m.rangereftofunc, cref)       
    
    delete_rows(m.inner, Vector{Int32}([r]))
end

function MOI.delete!(m::ClpMOIModel, cref::CR{MOI.SingleVariable, MOI.GreaterThan{Float64}})            
    lowerbounds = replaceInf(get_col_lower(m.inner))
    lowerbounds[m.lbreftovar[cref]] = -Inf    
    chg_column_lower(m.inner, lowerbounds)
    m.nlowerbounds -= 1
end

function MOI.delete!(m::ClpMOIModel, cref::CR{MOI.SingleVariable, MOI.LessThan{Float64}})            
    upperbounds = replaceInf(get_col_lower(m.inner))
    upperbounds[m.ubreftovar[cref]] = Inf    
    chg_column_upper(m.inner, upperbounds)
    m.nupperbounds -= 1
end

function MOI.delete!(m::ClpMOIModel, cref::CR{MOI.SingleVariable, MOI.EqualTo{Float64}})            
    upperbounds = replaceInf(get_col_lower(m.inner))
    upperbounds[m.ubreftovar[cref]] = Inf
    chg_column_upper(m.inner, upperbounds)
    lowerbounds = replaceInf(get_col_lower(m.inner))
    lowerbounds[m.lbreftovar[cref]] = -Inf    
    chg_column_lower(m.inner, lowerbounds)
    m.nassigns -= 1
end

function MOI.get(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.SingleVariable, MOI.LessThan{Float64}})
    return m.nupperbounds
end

function MOI.get(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.SingleVariable, MOI.GreaterThan{Float64}})
    return m.nlowerbounds
end

function MOI.get(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.SingleVariable, MOI.EqualTo{Float64}})
    return m.nassigns
end

function MOI.get(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}})
    return m.nrows_lessthan
end

function MOI.get(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.GreaterThan{Float64}})
    return m.nrows_greaterthan
end

function MOI.get(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.Interval{Float64}})
    return m.nrows_interval
end

function MOI.get(m::ClpMOIModel,
        ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}})
    return m.nrows_equalto
end

# Optimization and additional attributes
function MOI.optimize!(m::ClpMOIModel)
    initial_solve_with_options(m.inner,m.solveroptions)
    if is_proven_dual_infeasible(m.inner)
        # m.solution = unbounded_ray(m.inner) # .* get_obj_sense(m.inner)
    else
        m.solution = primal_column_solution(m.inner)
        m.rowsolution = primal_row_solution(m.inner)    
    end
    
    if is_proven_primal_infeasible(m.inner)
        # m.dualrowsolution = infeasibility_ray(m.inner) # .* get_obj_sense(m.inner)
    else
        m.dualrowsolution = dual_row_solution(m.inner) .* get_obj_sense(m.inner)
        m.dualcolsolution = dual_column_solution(m.inner) .* get_obj_sense(m.inner)
    end
    
    # @show MOI.get(m, MOI.ResultCount())
    # @show primal_feasible(m.inner)
    # @show dual_feasible(m.inner)
    # @show is_proven_primal_infeasible(m.inner)
    # @show is_proven_dual_infeasible(m.inner)
end

function MOI.get(m::ClpMOIModel, ::MOI.TerminationStatus)
    s = ClpCInterface.status(m.inner)
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
        error("Internal library error")
    end
end

function MOI.get(m::ClpMOIModel, ::MOI.PrimalStatus)
    # if is_proven_dual_infeasible(m.inner)
    #     return MOI.InfeasibilityCertificate        
    # else
    if primal_feasible(m.inner)
        return MOI.FeasiblePoint
    else    
        return MOI.UnknownResultStatus
    end
end

function MOI.get(m::ClpMOIModel, ::MOI.DualStatus)    
    # if is_proven_primal_infeasible(m.inner)
    #     return MOI.InfeasibilityCertificate       
    # else
    if dual_feasible(m.inner)
        return MOI.FeasiblePoint     
    else
        return MOI.UnknownResultStatus
    end
end

MOI.get(m::ClpMOIModel, ::MOI.ObjectiveValue) = objective_value(m.inner)

function MOI.get(m::ClpMOIModel, ::MOI.VariablePrimal, vref::VR)
    return m.solution[col(m,vref)+1]
end

function MOI.get(m::ClpMOIModel, ::MOI.VariablePrimal, vec::Vector{VR})
    return [m.solution[col(m,vref)+1] for vref in vec]
end

function MOI.get(m::ClpMOIModel, ::MOI.ConstraintPrimal, cref::CR)
    return m.rowsolution[row(m,cref)+1]
end

function MOI.get(m::ClpMOIModel, ::MOI.ConstraintPrimal, vec::Vector{CR})
    return [m.rowsolution[row(m,cref)+1] for cref in vec]
end

function MOI.get(m::ClpMOIModel, ::MOI.ConstraintDual, cref::CR{MOI.ScalarAffineFunction{Float64},S}
        where S <: Union{MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64}, MOI.EqualTo{Float64}})     
    return m.dualrowsolution[row(m,cref)+1]
end

function MOI.get(m::ClpMOIModel, ::MOI.ConstraintDual,
        cref::CR{MOI.SingleVariable, MOI.GreaterThan{Float64}})
    var = m.lbreftovar[cref]
    if m.solution[var] == get_col_lower(m.inner)[var] #this wouold be more efficient to get from instance
        return m.dualcolsolution[var]
    else
        return 0
    end
end

function MOI.get(m::ClpMOIModel, ::MOI.ConstraintDual,
        cref::CR{MOI.SingleVariable, MOI.LessThan{Float64}})
    var = m.ubreftovar[cref]
    if m.solution[var] == get_col_upper(m.inner)[var] #this wouold be more efficient to get from instance
        return m.dualcolsolution[var]
    else
        return 0
    end
end

function MOI.get(m::ClpMOIModel, ::MOI.ConstraintDual,
        cref::CR{MOI.SingleVariable, MOI.EqualTo{Float64}})
    var = m.assignreftovar[cref]    
    return m.dualcolsolution[var]    
end

function MOI.get(m::ClpMOIModel, ::MOI.ConstraintDual, vec::Vector{CR})
    return [MOI.get(m, MOI.ConstraintDual(), cref) for cref in vec]
end

MOI.canget(m::ClpMOIModel, ::Union{MOI.TerminationStatus,
                                    MOI.PrimalStatus,
                                    MOI.DualStatus,
                                    MOI.NumberOfVariables,
                                    MOI.NumberOfConstraints,
                                    MOI.ObjectiveSense,
                                    MOI.ObjectiveValue
                                    }) = true

MOI.canget(m::ClpMOIModel, ::MOI.PrimalStatus) = primal_feasible(m.inner) && !is_proven_primal_infeasible(m.inner)

MOI.canget(m::ClpMOIModel, ::MOI.DualStatus) = dual_feasible(m.inner) && !is_proven_dual_infeasible(m.inner)

MOI.canget(m::ClpMOIModel, ::Union{MOI.ConstraintPrimal, MOI.VariablePrimal, MOI.ConstraintDual},
        ::Union{VR, CR}) = true

MOI.canget(m::ClpMOIModel, ::Union{MOI.ConstraintPrimal, MOI.VariablePrimal},
        ::Union{Vector{VR}, Vector{CR}}) = true

MOI.canget(m::ClpMOIModel, ::MOI.ListOfConstraints) = false
MOI.canget(m::ClpMOIModel, ::Union{MOI.ConstraintFunction,
                                                 MOI.ConstraintSet}, ref::Union{VR, CR}) = false
MOI.canget(m::ClpMOIModel, ::MOI.ObjectiveFunction) = false
# MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintFunction, cr::CR)
# MOI.getattribute(m::ClpMOIModel, ::MOI.ConstraintSet, cr::CR)

MOI.canget(m::ClpMOIModel, ::MOI.ResultCount) = true

function MOI.get(m::ClpMOIModel, ::MOI.ResultCount) 
    if ((primal_feasible(m.inner) && !is_proven_dual_infeasible(m.inner)) 
            || (dual_feasible(m.inner) && !is_proven_primal_infeasible(m.inner)))
        return 1
    else
        return 0
    end    
end

MOI.get(solver::ClpMOISolver, ::Union{MOI.SupportsAddConstraintAfterSolve, MOI.SupportsAddVariableAfterSolve}) = true

MOI.get(solver::ClpMOISolver, ::MOI.SupportsDeleteConstraint) = true

MOI.get(solver::ClpMOISolver, ::MOI.SupportsDeleteVariable) = true

MOI.canmodifyconstraint(m::ClpMOIModel, ::CR, ::MOI.ScalarCoefficientChange{Float64}) = true

MOI.canmodifyconstraint(m::ClpMOIModel, ::CR, ::MOI.GreaterThan{Float64}) = true

MOI.canmodifyconstraint(m::ClpMOIModel, ::CR, ::MOI.LessThan{Float64}) = true

MOI.canmodifyconstraint(m::ClpMOIModel, ::CR, ::MOI.Interval{Float64}) = true

MOI.canmodifyobjective(m::ClpMOIModel, change::MOI.ScalarCoefficientChange{Float64}) = true

MOI.canmodifyobjective(m::ClpMOIModel, change::MOI.ScalarConstantChange{Float64}) = true

MOI.candelete(m::ClpMOIModel, ::CR{MOI.SingleVariable, MOI.GreaterThan{Float64}}) = true

MOI.candelete(m::ClpMOIModel, ::CR{MOI.SingleVariable, MOI.LessThan{Float64}}) = true

MOI.candelete(m::ClpMOIModel, ::CR{MOI.SingleVariable, MOI.EqualTo{Float64}}) = true

MOI.candelete(m::ClpMOIModel, ::CR{MOI.ScalarAffineFunction{Float64},S}) where S <: IGL = true

MOI.candelete(m::ClpMOIModel, ::VR) = true

end
