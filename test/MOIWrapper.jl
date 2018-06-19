using Base.Test, MathOptInterface, MathOptInterface.Test

const MOI  = MathOptInterface
const MOIT = MathOptInterface.Test

@testset "Linear tests" begin
    linconfig = MOIT.TestConfig(modify_lhs = false)
    @testset "Default Solver"  begin
        solver = ClpOptimizer(LogLevel = 0)
        MOIT.contlineartest(solver, linconfig, ["linear10","linear12","linear8a","linear8b","linear8c"])
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
    @testset "orderedindicestest" begin
        MOIT.orderedindicestest(solver)
    end
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

    MOIT.basic_constraint_tests(solver, config, exclude = [
        (MOI.SingleVariable, MOI.Integer),
        (MOI.SingleVariable, MOI.ZeroOne)        
    ])

    MOIT.unittest(solver, config, [
        "solve_affine_interval", 
        "solve_qcp_edge_cases",  
        "solve_qp_edge_cases"
    ])  
end