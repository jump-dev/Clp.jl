include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))

linprogtest(LPSolver(:Clp))
