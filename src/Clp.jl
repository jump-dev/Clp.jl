###
### COIN-OR Clp API Wrapper
###

isdefined(Base, :__precompile__) && __precompile__()

module Clp


include("ClpCInterface.jl")
include("ClpSolverInterface.jl")

using Clp.ClpMathProgSolverInterface
export ClpSolver

end # module
