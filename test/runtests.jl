using Clp

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(ClpSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
linprogsolvertest(ClpSolver())
