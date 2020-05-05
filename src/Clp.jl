###
### COIN-OR Clp API Wrapper
###

module Clp

include("C_interface.jl")
include("MPB_wrapper.jl")
include("MOI_wrapper.jl")

using Clp.ClpMathProgSolverInterface
export ClpSolver

end # module
