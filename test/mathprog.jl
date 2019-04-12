import MathProgBase
testdir = joinpath(dirname(pathof(MathProgBase)), "..", "test")

include(joinpath(testdir, "linprog.jl"))
linprogtest(ClpSolver())

include(joinpath(testdir, "linproginterface.jl"))
linprogsolvertest(ClpSolver())

include(joinpath(testdir, "conicinterface.jl"))
coniclineartest(ClpSolver())
