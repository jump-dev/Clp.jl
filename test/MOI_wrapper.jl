using Test
using MathOptInterface
import Clp

const MOI  = MathOptInterface

const OPTIMIZER = Clp.Optimizer()
MOI.set(OPTIMIZER, MOI.Silent(), true)

@testset "SolverName" begin
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "Clp"
end

@testset "supports_default_copy_to" begin
    @test !MOI.Utilities.supports_allocate_load(OPTIMIZER, false)
    @test !MOI.Utilities.supports_allocate_load(OPTIMIZER, true)
    @test !MOI.Utilities.supports_default_copy_to(OPTIMIZER, false)
    @test !MOI.Utilities.supports_default_copy_to(OPTIMIZER, true)
end

const CACHED = MOI.Utilities.CachingOptimizer(
    MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    OPTIMIZER
)

const BRIDGED = MOI.Bridges.full_bridge_optimizer(CACHED, Float64)

const CONFIG = MOI.Test.TestConfig(
    dual_objective_value = false,
)

@testset "basic_constraint_tests" begin
    MOI.Test.basic_constraint_tests(CACHED, CONFIG)
end

@testset "Unit Tests" begin
    MOI.Test.unittest(BRIDGED, CONFIG, [
        # Not supported by upstream.
        "number_threads",

        # Tests that require integer variables
        "solve_integer_edge_cases",
        "solve_zero_one_with_bounds_1",
        "solve_zero_one_with_bounds_2",
        "solve_zero_one_with_bounds_3",
        "solve_objbound_edge_cases",

        # Tests that require quadratic objective / constraints
        "solve_qcp_edge_cases",
        "solve_qp_edge_cases",
        "delete_soc_variables",
    ])
end

@testset "Linear tests" begin
    MOI.Test.contlineartest(BRIDGED, CONFIG, [
        # The linear12 test requires the InfeasibilityCertificate for variable
        # bounds. These are available through C++, but not via the C interface.
        "linear12",

        # MOI.VariablePrimalStart not supported.
        "partial_start"
    ])
end

@testset "ModelLike" begin
    @testset "nametest" begin
        MOI.Test.nametest(BRIDGED)
    end
    @testset "validtest" begin
        MOI.Test.validtest(BRIDGED)
    end
    @testset "emptytest" begin
        MOI.Test.emptytest(BRIDGED)
    end
end

@testset "Nonexistant unbounded ray" begin
    MOI.empty!(BRIDGED)
    x = MOI.add_variables(BRIDGED, 5)
    MOI.set(
        BRIDGED,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0)
    )
    MOI.set(BRIDGED, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.optimize!(BRIDGED)
    status = MOI.get(BRIDGED, MOI.TerminationStatus())
    @test status == MOI.DUAL_INFEASIBLE
end

@testset "RawParameter" begin
    model = Clp.Optimizer()
    MOI.set(model, MOI.RawParameter("LogLevel"), 1)
    @test MOI.get(model, MOI.RawParameter("LogLevel")) == 1
    @test MOI.get(model, MOI.RawParameter(:LogLevel)) == 1
    MOI.set(model, MOI.RawParameter(:LogLevel), 2)
    @test MOI.get(model, MOI.RawParameter("LogLevel")) == 2
    @test MOI.get(model, MOI.RawParameter(:LogLevel)) == 2

    MOI.set(model, MOI.RawParameter("SolveType"), 1)
    @test MOI.get(model, MOI.RawParameter("SolveType")) == 1
    @test MOI.get(model, MOI.RawParameter(:SolveType)) == 1
    MOI.set(model, MOI.RawParameter("SolveType"), 4)
    @test MOI.get(model, MOI.RawParameter("SolveType")) == 4
    @test MOI.get(model, MOI.RawParameter(:SolveType)) == 4

    MOI.set(model, MOI.RawParameter("PresolveType"), 1)
    @test MOI.get(model, MOI.RawParameter("PresolveType")) == 1
    @test MOI.get(model, MOI.RawParameter(:PresolveType)) == 1
    MOI.set(model, MOI.RawParameter("PresolveType"), 0)
    @test MOI.get(model, MOI.RawParameter("PresolveType")) == 0
    @test MOI.get(model, MOI.RawParameter(:PresolveType)) == 0
end

@testset "All parameters" begin
    model = Clp.Optimizer()
    param = MOI.RawParameter("NotAnOption")
    @test !MOI.supports(model, param)
    @test_throws MOI.UnsupportedAttribute(param) MOI.get(model, param)
    @test_throws MOI.UnsupportedAttribute(param) MOI.set(model, param, false)
    for key in Clp.SUPPORTED_PARAMETERS
        @test MOI.supports(model, MOI.RawParameter(key))
        value = MOI.get(model, MOI.RawParameter(key))
        MOI.set(model, MOI.RawParameter(key), value)
        @test MOI.get(model, MOI.RawParameter(key)) == value
    end
end

@testset "copy_to bug" begin
    model = MOI.Utilities.Model{Float64}()
    x = MOI.add_variable(model)
    con = [
        MOI.add_constraint(
            model,
            MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
            MOI.EqualTo(1.0),
        )
        for i = 1:2
    ]
    clp = Clp.Optimizer()
    index_map = MOI.copy_to(clp, model)
    @test index_map[con[1]] != index_map[con[2]]
end

@testset "Options after empty!" begin
    model = Clp.Optimizer()
    @test MOI.get(model, MOI.Silent()) == false
    MOI.set(model, MOI.Silent(), true)
    @test MOI.get(model, MOI.Silent()) == true
    MOI.empty!(model)
    @test MOI.get(model, MOI.Silent()) == true
end
