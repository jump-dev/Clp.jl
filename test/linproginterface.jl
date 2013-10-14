using Clp

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))

m = Clp.ClpSolverInterface.model()
linprogsolvertest(m)
