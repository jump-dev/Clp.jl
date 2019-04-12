using MathOptInterface

const MOI  = MathOptInterface
const MOIB = MathOptInterface.Bridges
const MOIT = MathOptInterface.Test

@testset "Unit Tests" begin
    config = MOIT.TestConfig()
    solver = Clp.Optimizer(LogLevel = 0)
    MOIT.basic_constraint_tests(solver, config)
    MOIT.unittest(solver, config, [
        "solve_qp_edge_cases",           # unsupported
        "solve_qcp_edge_cases",          # unsupported
        "solve_affine_interval",         # unsupported
        "solve_objbound_edge_cases",     # unsupported integer variables
        "solve_integer_edge_cases",      # unsupported integer variables
    ])
end

@testset "Linear tests" begin
    linconfig = MOIT.TestConfig(modify_lhs = false, basis = true)
    solver = Clp.Optimizer(LogLevel = 0)
    MOIT.contlineartest(solver, linconfig, [
        # linear1 test is disabled due to the following bug:
        # https://projects.coin-or.org/Clp/ticket/84
        "linear1",
        # linear10 test is tested below because it has interval sets.
        "linear10",
        "linear10b",
        # linear11 test is excluded as it fails on Linux for some reason.
        # It passes on Mac and Windows.
        "linear11",
        # linear12 test requires the InfeasibilityCertificate for variable
        # bounds. These are available through C++, but not via the C interface.
        "linear12",
        # partial_start requires VariablePrimalStart to be implemented by the
        # solver.
        "partial_start"
    ])

    @testset "Interval Bridge" begin
        MOIT.linear10test(MOIB.SplitInterval{Float64}(solver), linconfig)
        MOIT.linear10btest(MOIB.SplitInterval{Float64}(solver), linconfig)
    end
    @testset "Slack Bridge" begin
        MOIT.linear10test(MOIB.ScalarSlack{Float64}(solver), linconfig)
        MOIT.linear10btest(MOIB.ScalarSlack{Float64}(solver), linconfig)
    end
end

@testset "ModelLike tests" begin
    solver = Clp.Optimizer(LogLevel = 0)
    @testset "nametest" begin
        MOIT.nametest(solver)
    end
    @testset "validtest" begin
        MOIT.validtest(solver)
    end
    @testset "emptytest" begin
        MOIT.emptytest(solver)
    end
    @testset "copytest" begin
        solver2 = Clp.Optimizer(LogLevel = 0)
        MOIT.copytest(solver,solver2)
    end
    # Clp returns C_NULL when queried for the infeasibility ray in this case.
    @testset "Inexistant unbounded ray" begin
        o = Clp.Optimizer(LogLevel = 0)
        x = MOI.add_variables(o, 5)
        MOI.set(o, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
                   MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0))
        MOI.set(o, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        MOI.optimize!(o)
        status = MOI.get(o, MOI.TerminationStatus())
        @test status == MOI.DUAL_INFEASIBLE
    end
end
