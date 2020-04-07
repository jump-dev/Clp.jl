using Clp
using Compat
using Compat.Test
import Compat.Pkg

const MOI = Clp.MOI

# @testset "MathProgBase" begin
#     include("mathprog.jl")
# end
const C = Clp.ClpCInterface
function double_free_bug()
    model = C.@clp_ccall(newModel, Ptr{Cvoid}, ())
    options = C.@clpsolve_ccall(new, Ptr{Cvoid}, ())

    C.@clp_ccall(getNumRows, Cint, (Ptr{Cvoid},), model)
    C.@clp_ccall(getNumCols, Cint, (Ptr{Cvoid},), model)

    model = C.@clp_ccall(newModel, Ptr{Cvoid}, ())
    model = C.@clp_ccall(newModel, Ptr{Cvoid}, ())
    model = C.@clp_ccall(newModel, Ptr{Cvoid}, ())

    C.@clp_ccall(
        addColumns,
        Cvoid,
        (Ptr{Cvoid}, Int32, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
        model, 1, [-Inf], [Inf], [0.0], Int32[0, 0], Int32[], Float64[]
    )
    C.@clp_ccall(
        addColumns,
        Cvoid,
        (Ptr{Cvoid}, Int32, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
        model, 1, [-Inf], [Inf], [0.0], Int32[0, 0], Int32[], Float64[]
    )
    C.@clp_ccall(
        addColumns,
        Cvoid,
        (Ptr{Cvoid}, Int32, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
        model, 1, [-Inf], [Inf], [0.0], Int32[0, 0], Int32[], Float64[]
    )
    C.@clp_ccall(
        addColumns,
        Cvoid,
        (Ptr{Cvoid}, Int32, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
        model, 1, [-Inf], [Inf], [0.0], Int32[0, 0], Int32[], Float64[]
    )
    C.@clp_ccall(
        addColumns,
        Cvoid,
        (Ptr{Cvoid}, Int32, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}),
        model, 1, [-Inf], [Inf], [0.0], Int32[0, 0], Int32[], Float64[]
    )
    C.@clp_ccall(setObjSense, Cvoid, (Ptr{Cvoid}, Float64), model, -1)
    n = C.@clp_ccall(getNumCols, Cint, (Ptr{Cvoid},), model)
    C.@clp_ccall(chgObjCoefficients, Cvoid, (Ptr{Cvoid}, Ptr{Float64}), model, fill(1.0, n))
    C.@clp_ccall(setObjectiveOffset, Cvoid, (Ptr{Cvoid}, Float64), model, 0.0)
    C.@clp_ccall(initialSolveWithOptions, Cint, (Ptr{Cvoid},Ptr{Cvoid}), model, options)
    C.@clp_ccall(status, Cint, (Ptr{Cvoid},), model)
end

double_free_bug()
double_free_bug()
double_free_bug()
double_free_bug()

@testset "solve_unbounded_model (i)" begin
    model = Clp.MOIU.CachingOptimizer(
        Clp.MOIU.UniversalFallback(Clp.MOIU.Model{Float64}()),
        Clp.Optimizer()
    )

    config = Clp.MOI.Test.TestConfig()
    Clp.MOI.Test.solve_unbounded_model(model, config)
end

@testset "solve_unbounded_model (i)" begin
    model = Clp.MOIU.CachingOptimizer(
        Clp.MOIU.UniversalFallback(Clp.MOIU.Model{Float64}()),
        Clp.Optimizer()
    )
    MOI.empty!(model)
    x = MOI.add_variables(model, 5)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0)
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.DUAL_INFEASIBLE
end

# @testset "solve_unbounded_model (ii)" begin
#     model = Clp.MOIU.CachingOptimizer(
#         Clp.MOIU.UniversalFallback(Clp.MOIU.Model{Float64}()),
#         Clp.Optimizer()
#     )
#     config = Clp.MOI.Test.TestConfig()
#     Clp.MOI.Test.solve_unbounded_model(model, config)
# end

# @testset "solve_unbounded_model (iii)" begin
#     model = Clp.MOIU.CachingOptimizer(
#         Clp.MOIU.UniversalFallback(Clp.MOIU.Model{Float64}()),
#         Clp.Optimizer()
#     )
#     config = Clp.MOI.Test.TestConfig()
#     Clp.MOI.Test.solve_unbounded_model(model, config)
# end

# @testset "MathOptInterface" begin
#     include("MOIWrapper.jl")
# end
