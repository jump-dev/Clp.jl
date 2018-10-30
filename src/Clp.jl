###
### COIN-OR Clp API Wrapper
###

VERSION < v"0.7.0-beta2.199" && __precompile__()
module Clp

using Compat

include("ClpCInterface.jl")
include("ClpSolverInterface.jl")
include("MOIWrapper.jl")

using Clp.ClpMathProgSolverInterface
export ClpSolver

end # module
