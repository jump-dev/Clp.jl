###
### COIN-OR Clp API Wrapper
###

module Clp

using Compat
# This 'using' is required to suppress a warning about Clp not having Libdl in its
# dependencies (Libdl is used by BinaryProvider), e.g.: bicycle1885/CodecZlib.jl#26.
using Libdl

include("ClpCInterface.jl")
include("ClpSolverInterface.jl")
include("MOIWrapper.jl")

using Clp.ClpMathProgSolverInterface
export ClpSolver

end # module
