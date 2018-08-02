###
### COIN-OR Clp API Wrapper
###

__precompile__()

module Clp

import Pkg

include("ClpCInterface.jl")
include("ClpSolverInterface.jl")

using Clp.ClpMathProgSolverInterface

export ClpSolver

end # module
