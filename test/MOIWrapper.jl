using Test, MathOptInterface

const MOI  = MathOptInterface
const MOIB = MathOptInterface.Bridges
const MOIT = MathOptInterface.Test
const MOIU = MOI.Utilities

import Clp
const OPTIMIZER = Clp.Optimizer()
MOI.set(OPTIMIZER, MOI.Silent(), true)

@testset "SolverName" begin
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "Clp"
end

@testset "supports_default_copy_to" begin
    @test !MOIU.supports_allocate_load(OPTIMIZER, false)
    @test !MOIU.supports_allocate_load(OPTIMIZER, true)
    @test !MOIU.supports_default_copy_to(OPTIMIZER, false)
    @test !MOIU.supports_default_copy_to(OPTIMIZER, true)
end

const CACHE = MOIU.UniversalFallback(MOIU.Model{Float64}())
const CACHED = MOIU.CachingOptimizer(CACHE, OPTIMIZER)
MOI.set(CACHED, MOI.Silent(), true)

const BRIDGED = MOI.Bridges.full_bridge_optimizer(CACHED, Float64)
MOI.set(BRIDGED, MOI.Silent(), true)

const CONFIG = MOIT.TestConfig(
    dual_objective_value=false
)

@testset "basic_constraint_tests" begin
    MOIT.basic_constraint_tests(CACHED, CONFIG)
end

@testset "Unit Tests" begin
    MOIT.unittest(BRIDGED, CONFIG, [
        # Unsupported attributes
        "number_threads",  # not supported by Clp
        "time_limit_sec",  # Weird behaviour of Clp
        "solve_time",  # not supported by Clp
        # Tests that require integer variables
        "solve_integer_edge_cases",
        "solve_zero_one_with_bounds_1",
        "solve_zero_one_with_bounds_2",
        "solve_zero_one_with_bounds_3",
        "solve_objbound_edge_cases",
        # Tests that require quadratic objective / constraints
        "solve_qcp_edge_cases",
        "solve_qp_edge_cases",
        # Tests that require SOC
        "delete_soc_variables",
    ])
end


@testset "Linear tests" begin
    MOIT.contlineartest(BRIDGED, CONFIG, [
        # Queries an infeasibility certificate
        "linear8a",
        # linear12 test requires the InfeasibilityCertificate for variable
        # bounds. These are available through C++, but not via the C interface.
        "linear12",
        # partial_start requires VariablePrimalStart to be implemented by the
        # solver.
        "partial_start"
    ])
end

@testset "ModelLike" begin
    # solver = Clp.Optimizer(LogLevel = 0)
    @testset "nametest" begin
        MOIT.nametest(BRIDGED)
    end
    @testset "validtest" begin
        MOIT.validtest(BRIDGED)
    end
    @testset "emptytest" begin
        MOIT.emptytest(BRIDGED)
    end
    # @testset "copytest" begin
    #     solver2 = Clp.Optimizer(LogLevel = 0)
    #     MOIT.copytest(solver,solver2)
    # end
    # Clp returns C_NULL when queried for the infeasibility ray in this case.
    @testset "Inexistant unbounded ray" begin
        # o = Clp.Optimizer(LogLevel = 0)
        MOI.empty!(BRIDGED)

        x = MOI.add_variables(BRIDGED, 5)
        MOI.set(BRIDGED, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                   MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0))
        MOI.set(BRIDGED, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        MOI.optimize!(BRIDGED)
        status = MOI.get(BRIDGED, MOI.TerminationStatus())
        @test status == MOI.DUAL_INFEASIBLE
    end
end