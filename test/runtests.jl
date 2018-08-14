using Clp, Base.Test

@testset "MathProgBase" begin
    include("mathprog.jl")
end

@testset "MathOptInterface" begin
    include("MOIWrapper.jl")
end
