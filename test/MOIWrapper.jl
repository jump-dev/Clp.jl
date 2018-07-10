using Base.Test, MathOptInterface, MathOptInterface.Test

const MOI  = MathOptInterface
const MOIT = MathOptInterface.Test

@testset "Linear tests" begin
    linconfig = MOIT.TestConfig(modify_lhs = false)
    @testset "Default Solver"  begin
        solver = ClpOptimizer(LogLevel = 0)
        # linear1 test is disabled due to the following bug.
        # https://projects.coin-or.org/Clp/ticket/84
        MOIT.contlineartest(solver, linconfig, ["linear1", "linear10","linear11",
            "linear12","linear8a","linear8b","linear8c"])
    end
end

@testset "ModelLike tests" begin
    solver = ClpOptimizer(LogLevel = 0)
    MOIT.nametest(solver)
    @testset "validtest" begin
        MOIT.validtest(solver)
    end
    @testset "emptytest" begin
        MOIT.emptytest(solver)
    end
    # TODO uncomment after adding MOIT.TestConfig.multiple_bounds
    # @testset "orderedindicestest" begin
    #     MOIT.orderedindicestest(solver)
    # end
    @testset "canaddconstrainttest" begin
        MOIT.canaddconstrainttest(solver, Float64, Complex{Float64})
    end
    @testset "copytest" begin
        solver2 = ClpOptimizer(LogLevel = 0)
        MOIT.copytest(solver,solver2)
    end
end

@testset "Unit Tests" begin
    config = MOIT.TestConfig()
    solver = ClpOptimizer(LogLevel = 0)
    # TODO uncomment after adding MOIT.TestConfig.multiple_bounds
    # MOIT.basic_constraint_tests(solver, config)
    MOIT.unittest(solver, config, [
        "solve_qcp_edge_cases",  
        "solve_qp_edge_cases"
    ])  
end