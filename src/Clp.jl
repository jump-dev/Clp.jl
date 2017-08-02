###
### COIN-OR Clp API Wrapper
###

__precompile__()

module Clp


include("ClpCInterface.jl")
include("ClpSolverInterface.jl")

using Clp.ClpMathProgSolverInterface
using Clp.ClpMOISolverInterface

export ClpSolver
export ClpMOISolver

end # module
