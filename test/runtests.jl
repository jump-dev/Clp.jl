using Clp

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(ClpSolver(LogLevel=1))

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
linprogsolvertest(ClpSolver())
