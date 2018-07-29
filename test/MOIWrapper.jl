using MathOptInterface

const MOI  = MathOptInterface
const MOIB = MathOptInterface.Bridges
const MOIT = MathOptInterface.Test

@testset "Unit Tests" begin
    config = MOIT.TestConfig()
    solver = ClpOptimizer(LogLevel = 0)
    # MOIT.basic_constraint_tests(solver, config)
    MOIT.unittest(solver, config, [
        "solve_qcp_edge_cases",
        "solve_qp_edge_cases"
    ])
end

@testset "Linear tests" begin
    linconfig = MOIT.TestConfig(modify_lhs = false)
    solver = ClpOptimizer(LogLevel = 0)
    MOIT.contlineartest(solver, linconfig, [
        # linear1 test is disabled due to the following bug.
        # https://projects.coin-or.org/Clp/ticket/84
        "linear1",
        # linear10 test is tested below because it has interval sets
        "linear10",
        # linear12 test requires the InfeasibilityCertificate for variable
        # bounds. These are available through C++, but not via the C interface.
        "linear12"
    ])

    @testset "Interval Bridge" begin
        MOIT.linear10test(MOIB.SplitInterval{Float64}(solver), linconfig)
    end
end

@testset "ModelLike tests" begin
    solver = ClpOptimizer(LogLevel = 0)
    @testset "nametest" begin
        MOIT.nametest(solver)
    end
    @testset "validtest" begin
        MOIT.validtest(solver)
    end
    @testset "emptytest" begin
        MOIT.emptytest(solver)
    end
    # @testset "orderedindicestest" begin
        # MOIT.orderedindicestest(solver)
    # end
    @testset "canaddconstrainttest" begin
        MOIT.canaddconstrainttest(solver, Float64, Complex{Float64})
    end
    @testset "copytest" begin
        solver2 = ClpOptimizer(LogLevel = 0)
        MOIT.copytest(solver,solver2)
    end
end
