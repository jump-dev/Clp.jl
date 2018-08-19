using Clp
using Compat
using Compat.Test
using Compat.Pkg

@testset "MathProgBase" begin
    include("mathprog.jl")
end

@testset "MathOptInterface" begin
    include("MOIWrapper.jl")
end
