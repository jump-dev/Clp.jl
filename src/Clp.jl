###
### COIN-OR Clp API Wrapper
###


module Clp


include("ClpCInterface.jl")
include("ClpSolverInterface.jl")

using Clp.ClpMathProgSolverInterface
export ClpSolver

end # module
