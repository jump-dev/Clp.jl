using Clp

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))

linprogsolvertest(ClpLPSolver())
