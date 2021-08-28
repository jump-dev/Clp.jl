const Clp_Solve = Cvoid

"""
    Clp_Version()

Clp library version number as string. 
"""
function Clp_Version()
    ccall((:Clp_Version, libClp), Cstring, ())
end

"""
    Clp_VersionMajor()

Major number of Clp library version. 
"""
function Clp_VersionMajor()
    ccall((:Clp_VersionMajor, libClp), Cint, ())
end

"""
    Clp_VersionMinor()

Minor number of Clp library version. 
"""
function Clp_VersionMinor()
    ccall((:Clp_VersionMinor, libClp), Cint, ())
end

"""
    Clp_VersionRelease()

Release number of Clp library version. 
"""
function Clp_VersionRelease()
    ccall((:Clp_VersionRelease, libClp), Cint, ())
end

const Clp_Simplex = Cvoid

"""
    Clp_newModel()

Default constructor 
"""
function Clp_newModel()
    ccall((:Clp_newModel, libClp), Ptr{Clp_Simplex}, ())
end

"""
    Clp_deleteModel(model)

Destructor 
"""
function Clp_deleteModel(model)
    ccall((:Clp_deleteModel, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

# no prototype is found for this function at Clp_C_Interface.h:77:35, please use with caution
"""
    ClpSolve_new()

Default constructor 
"""
function ClpSolve_new()
    ccall((:ClpSolve_new, libClp), Ptr{Clp_Solve}, ())
end

"""
    ClpSolve_delete(solve)

Destructor 
"""
function ClpSolve_delete(solve)
    ccall((:ClpSolve_delete, libClp), Cvoid, (Ptr{Clp_Solve},), solve)
end

const CoinBigIndex = Cint

"""
    Clp_loadProblem(model, numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub)

Loads a problem (the constraints on the rows are given by lower and upper bounds). If a pointer is NULL then the following values are the default: <ul> <li> <code>colub</code>: all columns have upper bound infinity <li> <code>collb</code>: all columns have lower bound 0 <li> <code>rowub</code>: all rows have upper bound infinity <li> <code>rowlb</code>: all rows have lower bound -infinity <li> <code>obj</code>: all variables have 0 objective coefficient </ul>

Just like the other loadProblem() method except that the matrix is given in a standard column major ordered format (without gaps). 
"""
function Clp_loadProblem(model, numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub)
    ccall((:Clp_loadProblem, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cint, Ptr{CoinBigIndex}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), model, numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub)
end

function Clp_loadQuadraticObjective(model, numberColumns, start, column, element)
    ccall((:Clp_loadQuadraticObjective, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{CoinBigIndex}, Ptr{Cint}, Ptr{Cdouble}), model, numberColumns, start, column, element)
end

"""
    Clp_readMps(model, filename, keepNames, ignoreErrors)

Read an mps file from the given filename 
"""
function Clp_readMps(model, filename, keepNames, ignoreErrors)
    ccall((:Clp_readMps, libClp), Cint, (Ptr{Clp_Simplex}, Cstring, Cint, Cint), model, filename, keepNames, ignoreErrors)
end

"""
    Clp_writeMps(model, filename, formatType, numberAcross, objSense)

Write an mps file to the given filename 

Format type is 0 = normal, 1 = extra or 2 = hex. Number across is 1 or 2. Use objSense = -1D to flip the objective function around. 
"""
function Clp_writeMps(model, filename, formatType, numberAcross, objSense)
    ccall((:Clp_writeMps, libClp), Cint, (Ptr{Clp_Simplex}, Cstring, Cint, Cint, Cdouble), model, filename, formatType, numberAcross, objSense)
end

"""
    Clp_copyInIntegerInformation(model, information)

Copy in integer informations 
"""
function Clp_copyInIntegerInformation(model, information)
    ccall((:Clp_copyInIntegerInformation, libClp), Cvoid, (Ptr{Clp_Simplex}, Cstring), model, information)
end

"""
    Clp_deleteIntegerInformation(model)

Drop integer informations 
"""
function Clp_deleteIntegerInformation(model)
    ccall((:Clp_deleteIntegerInformation, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_resize(model, newNumberRows, newNumberColumns)

Resizes rim part of model 
"""
function Clp_resize(model, newNumberRows, newNumberColumns)
    ccall((:Clp_resize, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cint), model, newNumberRows, newNumberColumns)
end

"""
    Clp_deleteRows(model, number, which)

Deletes rows 
"""
function Clp_deleteRows(model, number, which)
    ccall((:Clp_deleteRows, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{Cint}), model, number, which)
end

"""
    Clp_addRows(model, number, rowLower, rowUpper, rowStarts, columns, elements)

Add rows 
"""
function Clp_addRows(model, number, rowLower, rowUpper, rowStarts, columns, elements)
    ccall((:Clp_addRows, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{CoinBigIndex}, Ptr{Cint}, Ptr{Cdouble}), model, number, rowLower, rowUpper, rowStarts, columns, elements)
end

"""
    Clp_deleteColumns(model, number, which)

Deletes columns 
"""
function Clp_deleteColumns(model, number, which)
    ccall((:Clp_deleteColumns, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{Cint}), model, number, which)
end

"""
    Clp_addColumns(model, number, columnLower, columnUpper, objective, columnStarts, rows, elements)

Add columns 
"""
function Clp_addColumns(model, number, columnLower, columnUpper, objective, columnStarts, rows, elements)
    ccall((:Clp_addColumns, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{CoinBigIndex}, Ptr{Cint}, Ptr{Cdouble}), model, number, columnLower, columnUpper, objective, columnStarts, rows, elements)
end

"""
    Clp_chgRowLower(model, rowLower)

Change row lower bounds 
"""
function Clp_chgRowLower(model, rowLower)
    ccall((:Clp_chgRowLower, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, rowLower)
end

"""
    Clp_chgRowUpper(model, rowUpper)

Change row upper bounds 
"""
function Clp_chgRowUpper(model, rowUpper)
    ccall((:Clp_chgRowUpper, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, rowUpper)
end

"""
    Clp_chgColumnLower(model, columnLower)

Change column lower bounds 
"""
function Clp_chgColumnLower(model, columnLower)
    ccall((:Clp_chgColumnLower, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, columnLower)
end

"""
    Clp_chgColumnUpper(model, columnUpper)

Change column upper bounds 
"""
function Clp_chgColumnUpper(model, columnUpper)
    ccall((:Clp_chgColumnUpper, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, columnUpper)
end

"""
    Clp_chgObjCoefficients(model, objIn)

Change objective coefficients 
"""
function Clp_chgObjCoefficients(model, objIn)
    ccall((:Clp_chgObjCoefficients, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, objIn)
end

function Clp_modifyCoefficient(model, row, column, newElement, keepZero)
    ccall((:Clp_modifyCoefficient, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cint, Cdouble, Bool), model, row, column, newElement, keepZero)
end

"""
    Clp_dropNames(model)

Drops names - makes lengthnames 0 and names empty 
"""
function Clp_dropNames(model)
    ccall((:Clp_dropNames, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_copyNames(model, rowNames, columnNames)

Copies in names 
"""
function Clp_copyNames(model, rowNames, columnNames)
    ccall((:Clp_copyNames, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cstring}, Ptr{Cstring}), model, rowNames, columnNames)
end

"""
    Clp_numberRows(model)

Number of rows 
"""
function Clp_numberRows(model)
    ccall((:Clp_numberRows, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_numberColumns(model)

Number of columns 
"""
function Clp_numberColumns(model)
    ccall((:Clp_numberColumns, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_primalTolerance(model)

Primal tolerance to use 
"""
function Clp_primalTolerance(model)
    ccall((:Clp_primalTolerance, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setPrimalTolerance(model, value)
    ccall((:Clp_setPrimalTolerance, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

"""
    Clp_dualTolerance(model)

Dual tolerance to use 
"""
function Clp_dualTolerance(model)
    ccall((:Clp_dualTolerance, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setDualTolerance(model, value)
    ccall((:Clp_setDualTolerance, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

"""
    Clp_dualObjectiveLimit(model)

Dual objective limit 
"""
function Clp_dualObjectiveLimit(model)
    ccall((:Clp_dualObjectiveLimit, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setDualObjectiveLimit(model, value)
    ccall((:Clp_setDualObjectiveLimit, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

"""
    Clp_objectiveOffset(model)

Objective offset 
"""
function Clp_objectiveOffset(model)
    ccall((:Clp_objectiveOffset, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setObjectiveOffset(model, value)
    ccall((:Clp_setObjectiveOffset, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

"""
    Clp_problemName(model, maxNumberCharacters, array)

Fills in array with problem name 
"""
function Clp_problemName(model, maxNumberCharacters, array)
    ccall((:Clp_problemName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, maxNumberCharacters, array)
end

function Clp_setProblemName(model, maxNumberCharacters, array)
    ccall((:Clp_setProblemName, libClp), Cint, (Ptr{Clp_Simplex}, Cint, Cstring), model, maxNumberCharacters, array)
end

"""
    Clp_numberIterations(model)

Number of iterations 
"""
function Clp_numberIterations(model)
    ccall((:Clp_numberIterations, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setNumberIterations(model, numberIterations)
    ccall((:Clp_setNumberIterations, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, numberIterations)
end

"""
    maximumIterations(model)

Maximum number of iterations 
"""
function maximumIterations(model)
    ccall((:maximumIterations, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setMaximumIterations(model, value)
    ccall((:Clp_setMaximumIterations, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, value)
end

"""
    Clp_maximumSeconds(model)

Maximum time in seconds (from when set called) 
"""
function Clp_maximumSeconds(model)
    ccall((:Clp_maximumSeconds, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setMaximumSeconds(model, value)
    ccall((:Clp_setMaximumSeconds, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

"""
    Clp_hitMaximumIterations(model)

Returns true if hit maximum iterations (or time) 
"""
function Clp_hitMaximumIterations(model)
    ccall((:Clp_hitMaximumIterations, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_status(model)

Status of problem: 0 - optimal 1 - primal infeasible 2 - dual infeasible 3 - stopped on iterations etc 4 - stopped due to errors
"""
function Clp_status(model)
    ccall((:Clp_status, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_setProblemStatus(model, problemStatus)

Set problem status 
"""
function Clp_setProblemStatus(model, problemStatus)
    ccall((:Clp_setProblemStatus, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, problemStatus)
end

"""
    Clp_secondaryStatus(model)

Secondary status of problem - may get extended 0 - none 1 - primal infeasible because dual limit reached 2 - scaled problem optimal - unscaled has primal infeasibilities 3 - scaled problem optimal - unscaled has dual infeasibilities 4 - scaled problem optimal - unscaled has both dual and primal infeasibilities
"""
function Clp_secondaryStatus(model)
    ccall((:Clp_secondaryStatus, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setSecondaryStatus(model, status)
    ccall((:Clp_setSecondaryStatus, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, status)
end

"""
    Clp_optimizationDirection(model)

Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore 
"""
function Clp_optimizationDirection(model)
    ccall((:Clp_optimizationDirection, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setOptimizationDirection(model, value)
    ccall((:Clp_setOptimizationDirection, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

"""
    Clp_primalRowSolution(model)

Primal row solution 
"""
function Clp_primalRowSolution(model)
    ccall((:Clp_primalRowSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_primalColumnSolution(model)

Primal column solution 
"""
function Clp_primalColumnSolution(model)
    ccall((:Clp_primalColumnSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_dualRowSolution(model)

Dual row solution 
"""
function Clp_dualRowSolution(model)
    ccall((:Clp_dualRowSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_dualColumnSolution(model)

Reduced costs 
"""
function Clp_dualColumnSolution(model)
    ccall((:Clp_dualColumnSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_rowLower(model)

Row lower 
"""
function Clp_rowLower(model)
    ccall((:Clp_rowLower, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_rowUpper(model)

Row upper 
"""
function Clp_rowUpper(model)
    ccall((:Clp_rowUpper, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_objective(model)

Objective 
"""
function Clp_objective(model)
    ccall((:Clp_objective, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_columnLower(model)

Column Lower 
"""
function Clp_columnLower(model)
    ccall((:Clp_columnLower, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_columnUpper(model)

Column Upper 
"""
function Clp_columnUpper(model)
    ccall((:Clp_columnUpper, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getNumElements(model)

Number of elements in matrix 
"""
function Clp_getNumElements(model)
    ccall((:Clp_getNumElements, libClp), CoinBigIndex, (Ptr{Clp_Simplex},), model)
end

function Clp_getVectorStarts(model)
    ccall((:Clp_getVectorStarts, libClp), Ptr{CoinBigIndex}, (Ptr{Clp_Simplex},), model)
end

function Clp_getIndices(model)
    ccall((:Clp_getIndices, libClp), Ptr{Cint}, (Ptr{Clp_Simplex},), model)
end

function Clp_getVectorLengths(model)
    ccall((:Clp_getVectorLengths, libClp), Ptr{Cint}, (Ptr{Clp_Simplex},), model)
end

function Clp_getElements(model)
    ccall((:Clp_getElements, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_objectiveValue(model)

Objective value 
"""
function Clp_objectiveValue(model)
    ccall((:Clp_objectiveValue, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_integerInformation(model)

Integer information 
"""
function Clp_integerInformation(model)
    ccall((:Clp_integerInformation, libClp), Cstring, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_infeasibilityRay(model)

Gives Infeasibility ray.

Use [`Clp_freeRay`](@ref) to free the returned array.

### Returns
infeasibility ray, or NULL returned if none/wrong.
"""
function Clp_infeasibilityRay(model)
    ccall((:Clp_infeasibilityRay, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_unboundedRay(model)

Gives ray in which the problem is unbounded.

Use [`Clp_freeRay`](@ref) to free the returned array.

### Returns
unbounded ray, or NULL returned if none/wrong.
"""
function Clp_unboundedRay(model)
    ccall((:Clp_unboundedRay, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_freeRay(model, ray)

Frees a infeasibility or unbounded ray. 
"""
function Clp_freeRay(model, ray)
    ccall((:Clp_freeRay, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, ray)
end

"""
    Clp_statusExists(model)

See if status array exists (partly for OsiClp) 
"""
function Clp_statusExists(model)
    ccall((:Clp_statusExists, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_statusArray(model)

Return address of status array (char[numberRows+numberColumns]) 
"""
function Clp_statusArray(model)
    ccall((:Clp_statusArray, libClp), Ptr{Cuchar}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_copyinStatus(model, statusArray)

Copy in status vector 
"""
function Clp_copyinStatus(model, statusArray)
    ccall((:Clp_copyinStatus, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cuchar}), model, statusArray)
end

function Clp_getColumnStatus(model, sequence)
    ccall((:Clp_getColumnStatus, libClp), Cint, (Ptr{Clp_Simplex}, Cint), model, sequence)
end

function Clp_getRowStatus(model, sequence)
    ccall((:Clp_getRowStatus, libClp), Cint, (Ptr{Clp_Simplex}, Cint), model, sequence)
end

function Clp_setColumnStatus(model, sequence, value)
    ccall((:Clp_setColumnStatus, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cint), model, sequence, value)
end

function Clp_setRowStatus(model, sequence, value)
    ccall((:Clp_setRowStatus, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cint), model, sequence, value)
end

"""
    Clp_setUserPointer(model, pointer)

User pointer for whatever reason 
"""
function Clp_setUserPointer(model, pointer)
    ccall((:Clp_setUserPointer, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cvoid}), model, pointer)
end

function Clp_getUserPointer(model)
    ccall((:Clp_getUserPointer, libClp), Ptr{Cvoid}, (Ptr{Clp_Simplex},), model)
end

# typedef void ( COINLINKAGE_CB * clp_callback ) ( Clp_Simplex * model , int msgno , int ndouble , const double * dvec , int nint , const int * ivec , int nchar , char * * cvec )
const clp_callback = Ptr{Cvoid}

"""
    Clp_registerCallBack(model, userCallBack)

Pass in Callback function. Message numbers up to 1000000 are Clp, Coin ones have 1000000 added 
"""
function Clp_registerCallBack(model, userCallBack)
    ccall((:Clp_registerCallBack, libClp), Cvoid, (Ptr{Clp_Simplex}, clp_callback), model, userCallBack)
end

"""
    Clp_clearCallBack(model)

Unset Callback function 
"""
function Clp_clearCallBack(model)
    ccall((:Clp_clearCallBack, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_setLogLevel(model, value)

Amount of print out: 0 - none 1 - just final 2 - just factorizations 3 - as 2 plus a bit more 4 - verbose above that 8,16,32 etc just for selective debug
"""
function Clp_setLogLevel(model, value)
    ccall((:Clp_setLogLevel, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, value)
end

function Clp_logLevel(model)
    ccall((:Clp_logLevel, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_lengthNames(model)

length of names (0 means no names0 
"""
function Clp_lengthNames(model)
    ccall((:Clp_lengthNames, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_rowName(model, iRow, name)

Fill in array (at least lengthNames+1 long) with a row name 
"""
function Clp_rowName(model, iRow, name)
    ccall((:Clp_rowName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, iRow, name)
end

"""
    Clp_columnName(model, iColumn, name)

Fill in array (at least lengthNames+1 long) with a column name 
"""
function Clp_columnName(model, iColumn, name)
    ccall((:Clp_columnName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, iColumn, name)
end

"""
    Clp_setRowName(model, iRow, name)

Set row name - Nice if they are short - 8 chars or less I think 
"""
function Clp_setRowName(model, iRow, name)
    ccall((:Clp_setRowName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, iRow, name)
end

"""
    Clp_setColumnName(model, iColumn, name)

Set column name - Nice if they are short - 8 chars or less I think 
"""
function Clp_setColumnName(model, iColumn, name)
    ccall((:Clp_setColumnName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, iColumn, name)
end

"""
    Clp_initialSolve(model)

General solve algorithm which can do presolve. See ClpSolve.hpp for options
"""
function Clp_initialSolve(model)
    ccall((:Clp_initialSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_initialSolveWithOptions(model, arg2)

Pass solve options. (Exception to direct analogue rule) 
"""
function Clp_initialSolveWithOptions(model, arg2)
    ccall((:Clp_initialSolveWithOptions, libClp), Cint, (Ptr{Clp_Simplex}, Ptr{Clp_Solve}), model, arg2)
end

"""
    Clp_initialDualSolve(model)

Dual initial solve 
"""
function Clp_initialDualSolve(model)
    ccall((:Clp_initialDualSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_initialPrimalSolve(model)

Primal initial solve 
"""
function Clp_initialPrimalSolve(model)
    ccall((:Clp_initialPrimalSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_initialBarrierSolve(model)

Barrier initial solve 
"""
function Clp_initialBarrierSolve(model)
    ccall((:Clp_initialBarrierSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_initialBarrierNoCrossSolve(model)

Barrier initial solve, no crossover 
"""
function Clp_initialBarrierNoCrossSolve(model)
    ccall((:Clp_initialBarrierNoCrossSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_dual(model, ifValuesPass)

Dual algorithm - see ClpSimplexDual.hpp for method 
"""
function Clp_dual(model, ifValuesPass)
    ccall((:Clp_dual, libClp), Cint, (Ptr{Clp_Simplex}, Cint), model, ifValuesPass)
end

"""
    Clp_primal(model, ifValuesPass)

Primal algorithm - see ClpSimplexPrimal.hpp for method 
"""
function Clp_primal(model, ifValuesPass)
    ccall((:Clp_primal, libClp), Cint, (Ptr{Clp_Simplex}, Cint), model, ifValuesPass)
end

"""
    Clp_idiot(model, tryhard)

Solve the problem with the idiot code 
"""
function Clp_idiot(model, tryhard)
    ccall((:Clp_idiot, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, tryhard)
end

"""
    Clp_scaling(model, mode)

Sets or unsets scaling, 0 -off, 1 equilibrium, 2 geometric, 3, auto, 4 dynamic(later) 
"""
function Clp_scaling(model, mode)
    ccall((:Clp_scaling, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, mode)
end

"""
    Clp_scalingFlag(model)

Gets scalingFlag 
"""
function Clp_scalingFlag(model)
    ccall((:Clp_scalingFlag, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_crash(model, gap, pivot)

Crash - at present just aimed at dual, returns -2 if dual preferred and crash basis created -1 if dual preferred and all slack basis preferred 0 if basis going in was not all slack 1 if primal preferred and all slack basis preferred 2 if primal preferred and crash basis created.

if gap between bounds <="gap" variables can be flipped

If "pivot" is 0 No pivoting (so will just be choice of algorithm) 1 Simple pivoting e.g. gub 2 Mini iterations
"""
function Clp_crash(model, gap, pivot)
    ccall((:Clp_crash, libClp), Cint, (Ptr{Clp_Simplex}, Cdouble, Cint), model, gap, pivot)
end

"""
    Clp_primalFeasible(model)

If problem is primal feasible 
"""
function Clp_primalFeasible(model)
    ccall((:Clp_primalFeasible, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_dualFeasible(model)

If problem is dual feasible 
"""
function Clp_dualFeasible(model)
    ccall((:Clp_dualFeasible, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_dualBound(model)

Dual bound 
"""
function Clp_dualBound(model)
    ccall((:Clp_dualBound, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setDualBound(model, value)
    ccall((:Clp_setDualBound, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

"""
    Clp_infeasibilityCost(model)

Infeasibility cost 
"""
function Clp_infeasibilityCost(model)
    ccall((:Clp_infeasibilityCost, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setInfeasibilityCost(model, value)
    ccall((:Clp_setInfeasibilityCost, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

"""
    Clp_perturbation(model)

Perturbation: 50 - switch on perturbation 100 - auto perturb if takes too long (1.0e-6 largest nonzero) 101 - we are perturbed 102 - don't try perturbing again default is 100 others are for playing
"""
function Clp_perturbation(model)
    ccall((:Clp_perturbation, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setPerturbation(model, value)
    ccall((:Clp_setPerturbation, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, value)
end

"""
    Clp_algorithm(model)

Current (or last) algorithm 
"""
function Clp_algorithm(model)
    ccall((:Clp_algorithm, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_setAlgorithm(model, value)

Set algorithm 
"""
function Clp_setAlgorithm(model, value)
    ccall((:Clp_setAlgorithm, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, value)
end

"""
    Clp_sumDualInfeasibilities(model)

Sum of dual infeasibilities 
"""
function Clp_sumDualInfeasibilities(model)
    ccall((:Clp_sumDualInfeasibilities, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_numberDualInfeasibilities(model)

Number of dual infeasibilities 
"""
function Clp_numberDualInfeasibilities(model)
    ccall((:Clp_numberDualInfeasibilities, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_sumPrimalInfeasibilities(model)

Sum of primal infeasibilities 
"""
function Clp_sumPrimalInfeasibilities(model)
    ccall((:Clp_sumPrimalInfeasibilities, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_numberPrimalInfeasibilities(model)

Number of primal infeasibilities 
"""
function Clp_numberPrimalInfeasibilities(model)
    ccall((:Clp_numberPrimalInfeasibilities, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_saveModel(model, fileName)

Save model to file, returns 0 if success. This is designed for use outside algorithms so does not save iterating arrays etc. It does not save any messaging information. Does not save scaling values. It does not know about all types of virtual functions.
"""
function Clp_saveModel(model, fileName)
    ccall((:Clp_saveModel, libClp), Cint, (Ptr{Clp_Simplex}, Cstring), model, fileName)
end

"""
    Clp_restoreModel(model, fileName)

Restore model from file, returns 0 if success, deletes current model 
"""
function Clp_restoreModel(model, fileName)
    ccall((:Clp_restoreModel, libClp), Cint, (Ptr{Clp_Simplex}, Cstring), model, fileName)
end

"""
    Clp_checkSolution(model)

Just check solution (for external use) - sets sum of infeasibilities etc 
"""
function Clp_checkSolution(model)
    ccall((:Clp_checkSolution, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getNumRows(model)

Number of rows 
"""
function Clp_getNumRows(model)
    ccall((:Clp_getNumRows, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getNumCols(model)

Number of columns 
"""
function Clp_getNumCols(model)
    ccall((:Clp_getNumCols, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getIterationCount(model)

Number of iterations 
"""
function Clp_getIterationCount(model)
    ccall((:Clp_getIterationCount, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_isAbandoned(model)

Are there a numerical difficulties? 
"""
function Clp_isAbandoned(model)
    ccall((:Clp_isAbandoned, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_isProvenOptimal(model)

Is optimality proven? 
"""
function Clp_isProvenOptimal(model)
    ccall((:Clp_isProvenOptimal, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_isProvenPrimalInfeasible(model)

Is primal infeasiblity proven? 
"""
function Clp_isProvenPrimalInfeasible(model)
    ccall((:Clp_isProvenPrimalInfeasible, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_isProvenDualInfeasible(model)

Is dual infeasiblity proven? 
"""
function Clp_isProvenDualInfeasible(model)
    ccall((:Clp_isProvenDualInfeasible, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_isPrimalObjectiveLimitReached(model)

Is the given primal objective limit reached? 
"""
function Clp_isPrimalObjectiveLimitReached(model)
    ccall((:Clp_isPrimalObjectiveLimitReached, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_isDualObjectiveLimitReached(model)

Is the given dual objective limit reached? 
"""
function Clp_isDualObjectiveLimitReached(model)
    ccall((:Clp_isDualObjectiveLimitReached, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_isIterationLimitReached(model)

Iteration limit reached? 
"""
function Clp_isIterationLimitReached(model)
    ccall((:Clp_isIterationLimitReached, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getObjSense(model)

Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore 
"""
function Clp_getObjSense(model)
    ccall((:Clp_getObjSense, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_setObjSense(model, objsen)

Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore 
"""
function Clp_setObjSense(model, objsen)
    ccall((:Clp_setObjSense, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, objsen)
end

"""
    Clp_getRowActivity(model)

Primal row solution 
"""
function Clp_getRowActivity(model)
    ccall((:Clp_getRowActivity, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getColSolution(model)

Primal column solution 
"""
function Clp_getColSolution(model)
    ccall((:Clp_getColSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_setColSolution(model, input)
    ccall((:Clp_setColSolution, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, input)
end

"""
    Clp_getRowPrice(model)

Dual row solution 
"""
function Clp_getRowPrice(model)
    ccall((:Clp_getRowPrice, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getReducedCost(model)

Reduced costs 
"""
function Clp_getReducedCost(model)
    ccall((:Clp_getReducedCost, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getRowLower(model)

Row lower 
"""
function Clp_getRowLower(model)
    ccall((:Clp_getRowLower, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getRowUpper(model)

Row upper 
"""
function Clp_getRowUpper(model)
    ccall((:Clp_getRowUpper, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getObjCoefficients(model)

Objective 
"""
function Clp_getObjCoefficients(model)
    ccall((:Clp_getObjCoefficients, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getColLower(model)

Column Lower 
"""
function Clp_getColLower(model)
    ccall((:Clp_getColLower, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getColUpper(model)

Column Upper 
"""
function Clp_getColUpper(model)
    ccall((:Clp_getColUpper, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_getObjValue(model)

Objective value 
"""
function Clp_getObjValue(model)
    ccall((:Clp_getObjValue, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

"""
    Clp_printModel(model, prefix)

Print model for debugging purposes 
"""
function Clp_printModel(model, prefix)
    ccall((:Clp_printModel, libClp), Cvoid, (Ptr{Clp_Simplex}, Cstring), model, prefix)
end

function Clp_getSmallElementValue(model)
    ccall((:Clp_getSmallElementValue, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setSmallElementValue(model, value)
    ccall((:Clp_setSmallElementValue, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

function ClpSolve_setSpecialOption(arg1, which, value, extraInfo)
    ccall((:ClpSolve_setSpecialOption, libClp), Cvoid, (Ptr{Clp_Solve}, Cint, Cint, Cint), arg1, which, value, extraInfo)
end

function ClpSolve_getSpecialOption(arg1, which)
    ccall((:ClpSolve_getSpecialOption, libClp), Cint, (Ptr{Clp_Solve}, Cint), arg1, which)
end

"""
    ClpSolve_setSolveType(arg1, method, extraInfo)

method: (see ClpSolve::SolveType) 0 - dual simplex 1 - primal simplex 2 - primal or sprint 3 - barrier 4 - barrier no crossover 5 - automatic 6 - not implemented -- pass extraInfo == -1 for default behavior 
"""
function ClpSolve_setSolveType(arg1, method, extraInfo)
    ccall((:ClpSolve_setSolveType, libClp), Cvoid, (Ptr{Clp_Solve}, Cint, Cint), arg1, method, extraInfo)
end

function ClpSolve_getSolveType(arg1)
    ccall((:ClpSolve_getSolveType, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

"""
    ClpSolve_setPresolveType(arg1, amount, extraInfo)

amount: (see ClpSolve::PresolveType) 0 - presolve on 1 - presolve off 2 - presolve number 3 - presolve number cost -- pass extraInfo == -1 for default behavior 
"""
function ClpSolve_setPresolveType(arg1, amount, extraInfo)
    ccall((:ClpSolve_setPresolveType, libClp), Cvoid, (Ptr{Clp_Solve}, Cint, Cint), arg1, amount, extraInfo)
end

function ClpSolve_getPresolveType(arg1)
    ccall((:ClpSolve_getPresolveType, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_getPresolvePasses(arg1)
    ccall((:ClpSolve_getPresolvePasses, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_getExtraInfo(arg1, which)
    ccall((:ClpSolve_getExtraInfo, libClp), Cint, (Ptr{Clp_Solve}, Cint), arg1, which)
end

function ClpSolve_setInfeasibleReturn(arg1, trueFalse)
    ccall((:ClpSolve_setInfeasibleReturn, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, trueFalse)
end

function ClpSolve_infeasibleReturn(arg1)
    ccall((:ClpSolve_infeasibleReturn, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_doDual(arg1)
    ccall((:ClpSolve_doDual, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoDual(arg1, doDual)
    ccall((:ClpSolve_setDoDual, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doDual)
end

function ClpSolve_doSingleton(arg1)
    ccall((:ClpSolve_doSingleton, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoSingleton(arg1, doSingleton)
    ccall((:ClpSolve_setDoSingleton, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doSingleton)
end

function ClpSolve_doDoubleton(arg1)
    ccall((:ClpSolve_doDoubleton, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoDoubleton(arg1, doDoubleton)
    ccall((:ClpSolve_setDoDoubleton, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doDoubleton)
end

function ClpSolve_doTripleton(arg1)
    ccall((:ClpSolve_doTripleton, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoTripleton(arg1, doTripleton)
    ccall((:ClpSolve_setDoTripleton, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doTripleton)
end

function ClpSolve_doTighten(arg1)
    ccall((:ClpSolve_doTighten, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoTighten(arg1, doTighten)
    ccall((:ClpSolve_setDoTighten, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doTighten)
end

function ClpSolve_doForcing(arg1)
    ccall((:ClpSolve_doForcing, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoForcing(arg1, doForcing)
    ccall((:ClpSolve_setDoForcing, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doForcing)
end

function ClpSolve_doImpliedFree(arg1)
    ccall((:ClpSolve_doImpliedFree, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoImpliedFree(arg1, doImpliedFree)
    ccall((:ClpSolve_setDoImpliedFree, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doImpliedFree)
end

function ClpSolve_doDupcol(arg1)
    ccall((:ClpSolve_doDupcol, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoDupcol(arg1, doDupcol)
    ccall((:ClpSolve_setDoDupcol, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doDupcol)
end

function ClpSolve_doDuprow(arg1)
    ccall((:ClpSolve_doDuprow, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoDuprow(arg1, doDuprow)
    ccall((:ClpSolve_setDoDuprow, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doDuprow)
end

function ClpSolve_doSingletonColumn(arg1)
    ccall((:ClpSolve_doSingletonColumn, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setDoSingletonColumn(arg1, doSingleton)
    ccall((:ClpSolve_setDoSingletonColumn, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, doSingleton)
end

function ClpSolve_presolveActions(arg1)
    ccall((:ClpSolve_presolveActions, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setPresolveActions(arg1, action)
    ccall((:ClpSolve_setPresolveActions, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, action)
end

function ClpSolve_substitution(arg1)
    ccall((:ClpSolve_substitution, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

function ClpSolve_setSubstitution(arg1, value)
    ccall((:ClpSolve_setSubstitution, libClp), Cvoid, (Ptr{Clp_Solve}, Cint), arg1, value)
end

const Sbb_Model = Cvoid

const Cbc_Model = Cvoid

# typedef void ( COINLINKAGE_CB * sbb_callback ) ( Sbb_Model * model , int msgno , int ndouble , const double * dvec , int nint , const int * ivec , int nchar , char * * cvec )
"""
typedef for user call back. The cvec are constructed so don't need to be const
"""
const sbb_callback = Ptr{Cvoid}

# typedef void ( COINLINKAGE_CB * cbc_callback ) ( Cbc_Model * model , int msgno , int ndouble , const double * dvec , int nint , const int * ivec , int nchar , char * * cvec )
const cbc_callback = Ptr{Cvoid}

# typedef void ( COINLINKAGE_CB * cbc_cut_callback ) ( void * osiSolver , void * osiCuts , void * appdata )
"""
typedef for cbc cut callback osiSolver needs to be an OsiSolverInterface object, osiCuts is an OsiCuts object and appdata is a pointer that will be passed to the cut  generation, you can use it to point to a data structure with information about the original problem,  for instance
"""
const cbc_cut_callback = Ptr{Cvoid}

