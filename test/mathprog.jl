if VERSION < v"0.7"
    testdir = joinpath(Pkg.dir("MathProgBase"), "test")
else
    import MathProgBase
    testdir = joinpath(dirname(pathof(MathProgBase)), "..", "test")
end
include(joinpath(testdir, "linprog.jl"))
linprogtest(ClpSolver())

include(joinpath(testdir, "linproginterface.jl"))
linprogsolvertest(ClpSolver())

include(joinpath(testdir, "conicinterface.jl"))
coniclineartest(ClpSolver())
