using Clp
using Compat
using Compat.Test
import Compat.Pkg

@testset "MathProgBase" begin
    include("mathprog.jl")
end

@testset "MathOptInterface" begin
    include("MOIWrapper.jl")
end
