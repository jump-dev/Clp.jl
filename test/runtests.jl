using Clp
using Compat
using Compat.Test
import Compat.Pkg

@testset "MathProgBase" begin
    include("MPB_wrapper.jl")
end

@testset "MathOptInterface" begin
    include("MOI_wrapper.jl")
end
