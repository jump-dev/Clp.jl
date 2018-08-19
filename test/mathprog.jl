# This is needed because both Compat.Pkg and MathProgBase export status :(.
# Good thing MOI's tests aren't loaded via include.
@static if VERSION < v"0.7-"
    using MathProgBase: status
end

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(ClpSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
linprogsolvertest(ClpSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
coniclineartest(ClpSolver())
