if get(ENV, "GITHUB_ACTIONS", "") == "true"
    import Pkg
    Pkg.add(Pkg.PackageSpec(name = "MathOptInterface", rev = "master"))
end

using Test

@testset "MathOptInterface" begin
    include("MOI_wrapper.jl")
end
