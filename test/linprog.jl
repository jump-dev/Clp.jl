using Clp

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))

linprogtest(ClpLPSolver())
