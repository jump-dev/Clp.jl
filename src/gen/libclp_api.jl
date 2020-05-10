# Julia wrapper for header: Coin_C_defines.h
# Automatically generated using Clang.jl

# Julia wrapper for header: Clp_C_Interface.h
# Automatically generated using Clang.jl


function Clp_Version()
    ccall((:Clp_Version, libClp), Cstring, ())
end

function Clp_VersionMajor()
    ccall((:Clp_VersionMajor, libClp), Cint, ())
end

function Clp_VersionMinor()
    ccall((:Clp_VersionMinor, libClp), Cint, ())
end

function Clp_VersionRelease()
    ccall((:Clp_VersionRelease, libClp), Cint, ())
end

function Clp_newModel()
    ccall((:Clp_newModel, libClp), Ptr{Clp_Simplex}, ())
end

function Clp_deleteModel(model)
    ccall((:Clp_deleteModel, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

function ClpSolve_new()
    ccall((:ClpSolve_new, libClp), Ptr{Clp_Solve}, ())
end

function ClpSolve_delete(solve)
    ccall((:ClpSolve_delete, libClp), Cvoid, (Ptr{Clp_Solve},), solve)
end

function Clp_loadProblem(model, numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub)
    ccall((:Clp_loadProblem, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cint, Ptr{CoinBigIndex}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), model, numcols, numrows, start, index, value, collb, colub, obj, rowlb, rowub)
end

function Clp_loadQuadraticObjective(model, numberColumns, start, column, element)
    ccall((:Clp_loadQuadraticObjective, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{CoinBigIndex}, Ptr{Cint}, Ptr{Cdouble}), model, numberColumns, start, column, element)
end

function Clp_readMps(model, filename, keepNames, ignoreErrors)
    ccall((:Clp_readMps, libClp), Cint, (Ptr{Clp_Simplex}, Cstring, Cint, Cint), model, filename, keepNames, ignoreErrors)
end

function Clp_writeMps(model, filename, formatType, numberAcross, objSense)
    ccall((:Clp_writeMps, libClp), Cint, (Ptr{Clp_Simplex}, Cstring, Cint, Cint, Cdouble), model, filename, formatType, numberAcross, objSense)
end

function Clp_copyInIntegerInformation(model, information)
    ccall((:Clp_copyInIntegerInformation, libClp), Cvoid, (Ptr{Clp_Simplex}, Cstring), model, information)
end

function Clp_deleteIntegerInformation(model)
    ccall((:Clp_deleteIntegerInformation, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

function Clp_resize(model, newNumberRows, newNumberColumns)
    ccall((:Clp_resize, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cint), model, newNumberRows, newNumberColumns)
end

function Clp_deleteRows(model, number, which)
    ccall((:Clp_deleteRows, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{Cint}), model, number, which)
end

function Clp_addRows(model, number, rowLower, rowUpper, rowStarts, columns, elements)
    ccall((:Clp_addRows, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{CoinBigIndex}, Ptr{Cint}, Ptr{Cdouble}), model, number, rowLower, rowUpper, rowStarts, columns, elements)
end

function Clp_deleteColumns(model, number, which)
    ccall((:Clp_deleteColumns, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{Cint}), model, number, which)
end

function Clp_addColumns(model, number, columnLower, columnUpper, objective, columnStarts, rows, elements)
    ccall((:Clp_addColumns, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{CoinBigIndex}, Ptr{Cint}, Ptr{Cdouble}), model, number, columnLower, columnUpper, objective, columnStarts, rows, elements)
end

function Clp_chgRowLower(model, rowLower)
    ccall((:Clp_chgRowLower, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, rowLower)
end

function Clp_chgRowUpper(model, rowUpper)
    ccall((:Clp_chgRowUpper, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, rowUpper)
end

function Clp_chgColumnLower(model, columnLower)
    ccall((:Clp_chgColumnLower, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, columnLower)
end

function Clp_chgColumnUpper(model, columnUpper)
    ccall((:Clp_chgColumnUpper, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, columnUpper)
end

function Clp_chgObjCoefficients(model, objIn)
    ccall((:Clp_chgObjCoefficients, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, objIn)
end

function Clp_modifyCoefficient(model, row, column, newElement, keepZero)
    ccall((:Clp_modifyCoefficient, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cint, Cdouble, Bool), model, row, column, newElement, keepZero)
end

function Clp_dropNames(model)
    ccall((:Clp_dropNames, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

function Clp_copyNames(model, rowNames, columnNames)
    ccall((:Clp_copyNames, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cstring}, Ptr{Cstring}), model, rowNames, columnNames)
end

function Clp_numberRows(model)
    ccall((:Clp_numberRows, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_numberColumns(model)
    ccall((:Clp_numberColumns, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_primalTolerance(model)
    ccall((:Clp_primalTolerance, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setPrimalTolerance(model, value)
    ccall((:Clp_setPrimalTolerance, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

function Clp_dualTolerance(model)
    ccall((:Clp_dualTolerance, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setDualTolerance(model, value)
    ccall((:Clp_setDualTolerance, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

function Clp_dualObjectiveLimit(model)
    ccall((:Clp_dualObjectiveLimit, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setDualObjectiveLimit(model, value)
    ccall((:Clp_setDualObjectiveLimit, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

function Clp_objectiveOffset(model)
    ccall((:Clp_objectiveOffset, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setObjectiveOffset(model, value)
    ccall((:Clp_setObjectiveOffset, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

function Clp_problemName(model, maxNumberCharacters, array)
    ccall((:Clp_problemName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, maxNumberCharacters, array)
end

function Clp_setProblemName(model, maxNumberCharacters, array)
    ccall((:Clp_setProblemName, libClp), Cint, (Ptr{Clp_Simplex}, Cint, Cstring), model, maxNumberCharacters, array)
end

function Clp_numberIterations(model)
    ccall((:Clp_numberIterations, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setNumberIterations(model, numberIterations)
    ccall((:Clp_setNumberIterations, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, numberIterations)
end

function maximumIterations(model)
    ccall((:maximumIterations, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setMaximumIterations(model, value)
    ccall((:Clp_setMaximumIterations, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, value)
end

function Clp_maximumSeconds(model)
    ccall((:Clp_maximumSeconds, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setMaximumSeconds(model, value)
    ccall((:Clp_setMaximumSeconds, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

function Clp_hitMaximumIterations(model)
    ccall((:Clp_hitMaximumIterations, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_status(model)
    ccall((:Clp_status, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setProblemStatus(model, problemStatus)
    ccall((:Clp_setProblemStatus, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, problemStatus)
end

function Clp_secondaryStatus(model)
    ccall((:Clp_secondaryStatus, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setSecondaryStatus(model, status)
    ccall((:Clp_setSecondaryStatus, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, status)
end

function Clp_optimizationDirection(model)
    ccall((:Clp_optimizationDirection, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setOptimizationDirection(model, value)
    ccall((:Clp_setOptimizationDirection, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

function Clp_primalRowSolution(model)
    ccall((:Clp_primalRowSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_primalColumnSolution(model)
    ccall((:Clp_primalColumnSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_dualRowSolution(model)
    ccall((:Clp_dualRowSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_dualColumnSolution(model)
    ccall((:Clp_dualColumnSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_rowLower(model)
    ccall((:Clp_rowLower, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_rowUpper(model)
    ccall((:Clp_rowUpper, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_objective(model)
    ccall((:Clp_objective, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_columnLower(model)
    ccall((:Clp_columnLower, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_columnUpper(model)
    ccall((:Clp_columnUpper, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

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

function Clp_objectiveValue(model)
    ccall((:Clp_objectiveValue, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_integerInformation(model)
    ccall((:Clp_integerInformation, libClp), Cstring, (Ptr{Clp_Simplex},), model)
end

function Clp_infeasibilityRay(model)
    ccall((:Clp_infeasibilityRay, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_unboundedRay(model)
    ccall((:Clp_unboundedRay, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_freeRay(model, ray)
    ccall((:Clp_freeRay, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, ray)
end

function Clp_statusExists(model)
    ccall((:Clp_statusExists, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_statusArray(model)
    ccall((:Clp_statusArray, libClp), Ptr{Cuchar}, (Ptr{Clp_Simplex},), model)
end

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

function Clp_setUserPointer(model, pointer)
    ccall((:Clp_setUserPointer, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cvoid}), model, pointer)
end

function Clp_getUserPointer(model)
    ccall((:Clp_getUserPointer, libClp), Ptr{Cvoid}, (Ptr{Clp_Simplex},), model)
end

function Clp_registerCallBack(model, userCallBack)
    ccall((:Clp_registerCallBack, libClp), Cvoid, (Ptr{Clp_Simplex}, clp_callback), model, userCallBack)
end

function Clp_clearCallBack(model)
    ccall((:Clp_clearCallBack, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

function Clp_setLogLevel(model, value)
    ccall((:Clp_setLogLevel, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, value)
end

function Clp_logLevel(model)
    ccall((:Clp_logLevel, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_lengthNames(model)
    ccall((:Clp_lengthNames, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_rowName(model, iRow, name)
    ccall((:Clp_rowName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, iRow, name)
end

function Clp_columnName(model, iColumn, name)
    ccall((:Clp_columnName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, iColumn, name)
end

function Clp_setRowName(model, iRow, name)
    ccall((:Clp_setRowName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, iRow, name)
end

function Clp_setColumnName(model, iColumn, name)
    ccall((:Clp_setColumnName, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint, Cstring), model, iColumn, name)
end

function Clp_initialSolve(model)
    ccall((:Clp_initialSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_initialSolveWithOptions(model, arg1)
    ccall((:Clp_initialSolveWithOptions, libClp), Cint, (Ptr{Clp_Simplex}, Ptr{Clp_Solve}), model, arg1)
end

function Clp_initialDualSolve(model)
    ccall((:Clp_initialDualSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_initialPrimalSolve(model)
    ccall((:Clp_initialPrimalSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_initialBarrierSolve(model)
    ccall((:Clp_initialBarrierSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_initialBarrierNoCrossSolve(model)
    ccall((:Clp_initialBarrierNoCrossSolve, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_dual(model, ifValuesPass)
    ccall((:Clp_dual, libClp), Cint, (Ptr{Clp_Simplex}, Cint), model, ifValuesPass)
end

function Clp_primal(model, ifValuesPass)
    ccall((:Clp_primal, libClp), Cint, (Ptr{Clp_Simplex}, Cint), model, ifValuesPass)
end

function Clp_idiot(model, tryhard)
    ccall((:Clp_idiot, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, tryhard)
end

function Clp_scaling(model, mode)
    ccall((:Clp_scaling, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, mode)
end

function Clp_scalingFlag(model)
    ccall((:Clp_scalingFlag, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_crash(model, gap, pivot)
    ccall((:Clp_crash, libClp), Cint, (Ptr{Clp_Simplex}, Cdouble, Cint), model, gap, pivot)
end

function Clp_primalFeasible(model)
    ccall((:Clp_primalFeasible, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_dualFeasible(model)
    ccall((:Clp_dualFeasible, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_dualBound(model)
    ccall((:Clp_dualBound, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setDualBound(model, value)
    ccall((:Clp_setDualBound, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

function Clp_infeasibilityCost(model)
    ccall((:Clp_infeasibilityCost, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setInfeasibilityCost(model, value)
    ccall((:Clp_setInfeasibilityCost, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, value)
end

function Clp_perturbation(model)
    ccall((:Clp_perturbation, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setPerturbation(model, value)
    ccall((:Clp_setPerturbation, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, value)
end

function Clp_algorithm(model)
    ccall((:Clp_algorithm, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_setAlgorithm(model, value)
    ccall((:Clp_setAlgorithm, libClp), Cvoid, (Ptr{Clp_Simplex}, Cint), model, value)
end

function Clp_sumDualInfeasibilities(model)
    ccall((:Clp_sumDualInfeasibilities, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_numberDualInfeasibilities(model)
    ccall((:Clp_numberDualInfeasibilities, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_sumPrimalInfeasibilities(model)
    ccall((:Clp_sumPrimalInfeasibilities, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_numberPrimalInfeasibilities(model)
    ccall((:Clp_numberPrimalInfeasibilities, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_saveModel(model, fileName)
    ccall((:Clp_saveModel, libClp), Cint, (Ptr{Clp_Simplex}, Cstring), model, fileName)
end

function Clp_restoreModel(model, fileName)
    ccall((:Clp_restoreModel, libClp), Cint, (Ptr{Clp_Simplex}, Cstring), model, fileName)
end

function Clp_checkSolution(model)
    ccall((:Clp_checkSolution, libClp), Cvoid, (Ptr{Clp_Simplex},), model)
end

function Clp_getNumRows(model)
    ccall((:Clp_getNumRows, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_getNumCols(model)
    ccall((:Clp_getNumCols, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_getIterationCount(model)
    ccall((:Clp_getIterationCount, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_isAbandoned(model)
    ccall((:Clp_isAbandoned, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_isProvenOptimal(model)
    ccall((:Clp_isProvenOptimal, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_isProvenPrimalInfeasible(model)
    ccall((:Clp_isProvenPrimalInfeasible, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_isProvenDualInfeasible(model)
    ccall((:Clp_isProvenDualInfeasible, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_isPrimalObjectiveLimitReached(model)
    ccall((:Clp_isPrimalObjectiveLimitReached, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_isDualObjectiveLimitReached(model)
    ccall((:Clp_isDualObjectiveLimitReached, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_isIterationLimitReached(model)
    ccall((:Clp_isIterationLimitReached, libClp), Cint, (Ptr{Clp_Simplex},), model)
end

function Clp_getObjSense(model)
    ccall((:Clp_getObjSense, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

function Clp_setObjSense(model, objsen)
    ccall((:Clp_setObjSense, libClp), Cvoid, (Ptr{Clp_Simplex}, Cdouble), model, objsen)
end

function Clp_getRowActivity(model)
    ccall((:Clp_getRowActivity, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_getColSolution(model)
    ccall((:Clp_getColSolution, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_setColSolution(model, input)
    ccall((:Clp_setColSolution, libClp), Cvoid, (Ptr{Clp_Simplex}, Ptr{Cdouble}), model, input)
end

function Clp_getRowPrice(model)
    ccall((:Clp_getRowPrice, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_getReducedCost(model)
    ccall((:Clp_getReducedCost, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_getRowLower(model)
    ccall((:Clp_getRowLower, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_getRowUpper(model)
    ccall((:Clp_getRowUpper, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_getObjCoefficients(model)
    ccall((:Clp_getObjCoefficients, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_getColLower(model)
    ccall((:Clp_getColLower, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_getColUpper(model)
    ccall((:Clp_getColUpper, libClp), Ptr{Cdouble}, (Ptr{Clp_Simplex},), model)
end

function Clp_getObjValue(model)
    ccall((:Clp_getObjValue, libClp), Cdouble, (Ptr{Clp_Simplex},), model)
end

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

function ClpSolve_setSolveType(arg1, method, extraInfo)
    ccall((:ClpSolve_setSolveType, libClp), Cvoid, (Ptr{Clp_Solve}, Cint, Cint), arg1, method, extraInfo)
end

function ClpSolve_getSolveType(arg1)
    ccall((:ClpSolve_getSolveType, libClp), Cint, (Ptr{Clp_Solve},), arg1)
end

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
