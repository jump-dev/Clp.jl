using Clp
using Compat
using Compat.Test
import Compat.Pkg

@testset "MathProgBase" begin
    include("mathprog.jl")
end

@testset "solve_unbounded_model (i)" begin
    model = Clp.MOIU.CachingOptimizer(
        Clp.MOIU.UniversalFallback(Clp.MOIU.Model{Float64}()),
        Clp.Optimizer()
    )
    config = Clp.MOI.Test.TestConfig()
    Clp.MOI.Test.solve_unbounded_model(model, config)
end

@testset "solve_unbounded_model (ii)" begin
    model = Clp.MOIU.CachingOptimizer(
        Clp.MOIU.UniversalFallback(Clp.MOIU.Model{Float64}()),
        Clp.Optimizer()
    )
    config = Clp.MOI.Test.TestConfig()
    Clp.MOI.Test.solve_unbounded_model(model, config)
end

@testset "solve_unbounded_model (iii)" begin
    model = Clp.MOIU.CachingOptimizer(
        Clp.MOIU.UniversalFallback(Clp.MOIU.Model{Float64}()),
        Clp.Optimizer()
    )
    config = Clp.MOI.Test.TestConfig()
    Clp.MOI.Test.solve_unbounded_model(model, config)
end

@testset "MathOptInterface" begin
    include("MOIWrapper.jl")
end
