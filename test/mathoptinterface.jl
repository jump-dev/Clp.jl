using MathOptInterface, Clp, Base.Test
const MOI = MathOptInterface

solver = ClpMOISolver()
ε=Base.rtoldefault(Float64)

function linear1test(solver::MOI.AbstractSolver, ε=Base.rtoldefault(Float64))
    @testset "Basic solve, query, resolve" begin
        # simple 2 variable, 1 constraint problem
        # min -x
        # st   x + y <= 1   (x + y - 1 ∈ Nonpositives)
        #       x, y >= 0   (x, y ∈ Nonnegatives)

        @test MOI.supportsproblem(solver, MOI.ScalarAffineFunction{Float64}, [(MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}),(MOI.SingleVariable,MOI.GreaterThan{Float64})])

        m = MOI.SolverInstance(solver)

        v = MOI.addvariables!(m, 2)
        @test MOI.getattribute(m, MOI.NumberOfVariables()) == 2

        cf = MOI.ScalarAffineFunction(v, [1.0,1.0], 0.0)
        c = MOI.addconstraint!(m, cf, MOI.LessThan(1.0))
        @test MOI.getattribute(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}}()) == 1

        vc1 = MOI.addconstraint!(m, MOI.SingleVariable(v[1]), MOI.GreaterThan(0.0))
        vc2 = MOI.addconstraint!(m, MOI.SingleVariable(v[2]), MOI.GreaterThan(0.0))
        @test MOI.getattribute(m, MOI.NumberOfConstraints{MOI.SingleVariable,MOI.GreaterThan{Float64}}()) == 2

        objf = MOI.ScalarAffineFunction(v, [-1.0,0.0], 0.0)
        MOI.setobjective!(m, MOI.MinSense, objf)

        @test MOI.getattribute(m, MOI.ObjectiveSense()) == MOI.MinSense

        if MOI.cangetattribute(m, MOI.ObjectiveFunction())
            @test objf ≈ MOI.getattribute(m, MOI.ObjectiveFunction())
        end

        if MOI.cangetattribute(m, MOI.ConstraintFunction(), c)
            @test cf ≈ MOI.getattribute(m, MOI.ConstraintFunction(), c)
        end

        if MOI.cangetattribute(m, MOI.ConstraintSet(), c)
            s = MOI.getattribute(m, MOI.ConstraintSet(), c)
            @test s == MOI.LessThan(1.0)
        end

        if MOI.cangetattribute(m, MOI.ConstraintSet(), vc1)
            s = MOI.getattribute(m, MOI.ConstraintSet(), vc1)
            @test s == MOI.GreaterThan(0.0)
        end
        if MOI.cangetattribute(m, MOI.ConstraintSet(), vc2)
            s = MOI.getattribute(m, MOI.ConstraintSet(), vc2)
            @test s == MOI.GreaterThan(0.0)
        end

        MOI.optimize!(m)

        @test MOI.cangetattribute(m, MOI.TerminationStatus())
        @test MOI.getattribute(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.cangetattribute(m, MOI.PrimalStatus())
        @test MOI.getattribute(m, MOI.PrimalStatus()) == MOI.FeasiblePoint

        @test MOI.cangetattribute(m, MOI.ObjectiveValue())
        @test MOI.getattribute(m, MOI.ObjectiveValue()) ≈ -1 atol=ε

        @test MOI.cangetattribute(m, MOI.VariablePrimal(), v)
        @test MOI.getattribute(m, MOI.VariablePrimal(), v) ≈ [1, 0] atol=ε

        @test MOI.cangetattribute(m, MOI.ConstraintPrimal(), c)
        @test MOI.getattribute(m, MOI.ConstraintPrimal(), c) ≈ 1 atol=ε

        if MOI.getattribute(solver, MOI.SupportsDuals())
            @test MOI.cangetattribute(m, MOI.DualStatus())
            @test MOI.getattribute(m, MOI.DualStatus()) == MOI.FeasiblePoint
            @test MOI.cangetattribute(m, MOI.ConstraintDual(), c)
            @test MOI.getattribute(m, MOI.ConstraintDual(), c) ≈ -1 atol=ε

            # reduced costs
            @test MOI.cangetattribute(m, MOI.ConstraintDual(), vc1)
            @test MOI.getattribute(m, MOI.ConstraintDual(), vc1) ≈ 0 atol=ε
            @test MOI.cangetattribute(m, MOI.ConstraintDual(), vc2)
            @test MOI.getattribute(m, MOI.ConstraintDual(), vc2) ≈ 1 atol=ε
        end

        # change objective to Max +x

        objf = MOI.ScalarAffineFunction(v, [1.0,0.0], 0.0)
        MOI.setobjective!(m, MOI.MaxSense, objf)

        if MOI.cangetattribute(m, MOI.ObjectiveFunction())
            @test objf ≈ MOI.getattribute(m, MOI.ObjectiveFunction())
        end

        @test MOI.getattribute(m, MOI.ObjectiveSense()) == MOI.MaxSense

        MOI.optimize!(m)

        @test MOI.cangetattribute(m, MOI.TerminationStatus())
        @test MOI.getattribute(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.cangetattribute(m, MOI.PrimalStatus())
        @test MOI.getattribute(m, MOI.PrimalStatus()) == MOI.FeasiblePoint

        @test MOI.cangetattribute(m, MOI.ObjectiveValue())
        @test MOI.getattribute(m, MOI.ObjectiveValue()) ≈ 1 atol=ε

        @test MOI.cangetattribute(m, MOI.VariablePrimal(), v)
        @test MOI.getattribute(m, MOI.VariablePrimal(), v) ≈ [1, 0] atol=ε

        if MOI.getattribute(solver, MOI.SupportsDuals())
            @test MOI.cangetattribute(m, MOI.DualStatus())
            @test MOI.getattribute(m, MOI.DualStatus()) == MOI.FeasiblePoint
            @test MOI.cangetattribute(m, MOI.ConstraintDual(), c)
            @test MOI.getattribute(m, MOI.ConstraintDual(), c) ≈ -1 atol=ε

            @test MOI.cangetattribute(m, MOI.ConstraintDual(), vc1)
            @test MOI.getattribute(m, MOI.ConstraintDual(), vc1) ≈ 0 atol=ε
            @test MOI.cangetattribute(m, MOI.ConstraintDual(), vc2)
            @test MOI.getattribute(m, MOI.ConstraintDual(), vc2) ≈ 1 atol=ε
        end

        # add new variable to get :
        # max x + 2z
        # s.t. x + y + z <= 1
        # x,y,z >= 0

        z = MOI.addvariable!(m)
        push!(v, z)
        @test v[3] == z

        vc3 = MOI.addconstraint!(m, MOI.SingleVariable(v[3]), MOI.GreaterThan(0.0))
        @test MOI.getattribute(m, MOI.NumberOfConstraints{MOI.SingleVariable,MOI.GreaterThan{Float64}}()) == 3

        # MOI.modifyconstraint!(m, c, MOI.ScalarCoefficientChange{Float64}(z, 1.0))
        # MOI.modifyobjective!(m, MOI.ScalarCoefficientChange{Float64}(z, 2.0))


    end
end

linear1test(solver, ε)
