using MathOptInterface, Clp, Base.Test
const MOI = MathOptInterface

solver = ClpMOISolver()
ε=Base.rtoldefault(Float64)


function linear1test(solver::MOI.AbstractSolver; atol=Base.rtoldefault(Float64), rtol=Base.rtoldefault(Float64))
    @testset "Basic solve, query, resolve" begin
        # simple 2 variable, 1 constraint problem
        # min -x
        # st   x + y <= 1   (x + y - 1 ∈ Nonpositives)
        #       x, y >= 0   (x, y ∈ Nonnegatives)

        @test MOI.supportsproblem(solver, MOI.ScalarAffineFunction{Float64}, [(MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}),(MOI.SingleVariable,MOI.GreaterThan{Float64})])

        @test MOI.get(solver, MOI.SupportsAddConstraintAfterSolve())
        @test MOI.get(solver, MOI.SupportsAddVariableAfterSolve())
        @test MOI.get(solver, MOI.SupportsDeleteConstraint())

        m = MOI.SolverInstance(solver)

        v = MOI.addvariables!(m, 2)
        @test MOI.get(m, MOI.NumberOfVariables()) == 2

        cf = MOI.ScalarAffineFunction(v, [1.0,1.0], 0.0)
        c = MOI.addconstraint!(m, cf, MOI.LessThan(1.0))
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}}()) == 1

        vc1 = MOI.addconstraint!(m, MOI.SingleVariable(v[1]), MOI.GreaterThan(0.0))
        # test fallback
        vc2 = MOI.addconstraint!(m, v[2], MOI.GreaterThan(0.0))
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.SingleVariable,MOI.GreaterThan{Float64}}()) == 2

        objf = MOI.ScalarAffineFunction(v, [-1.0,0.0], 0.0)
        MOI.set!(m, MOI.ObjectiveFunction(), objf)
        MOI.set!(m, MOI.ObjectiveSense(), MOI.MinSense)

        @test MOI.get(m, MOI.ObjectiveSense()) == MOI.MinSense

        if MOI.canget(m, MOI.ObjectiveFunction())
            @test objf ≈ MOI.get(m, MOI.ObjectiveFunction())
        end

        if MOI.canget(m, MOI.ConstraintFunction(), c)
            @test cf ≈ MOI.get(m, MOI.ConstraintFunction(), c)
        end

        if MOI.canget(m, MOI.ConstraintSet(), c)
            s = MOI.get(m, MOI.ConstraintSet(), c)
            @test s == MOI.LessThan(1.0)
        end

        if MOI.canget(m, MOI.ConstraintSet(), vc1)
            s = MOI.get(m, MOI.ConstraintSet(), vc1)
            @test s == MOI.GreaterThan(0.0)
        end
        if MOI.canget(m, MOI.ConstraintSet(), vc2)
            s = MOI.get(m, MOI.ConstraintSet(), vc2)
            @test s == MOI.GreaterThan(0.0)
        end

        MOI.optimize!(m)

        @test MOI.canget(m, MOI.TerminationStatus())
        @test MOI.get(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.canget(m, MOI.PrimalStatus())
        @test MOI.get(m, MOI.PrimalStatus()) == MOI.FeasiblePoint

        @test MOI.canget(m, MOI.ObjectiveValue())
        @test MOI.get(m, MOI.ObjectiveValue()) ≈ -1 atol=atol rtol=rtol

        @test MOI.canget(m, MOI.VariablePrimal(), v)
        @test MOI.get(m, MOI.VariablePrimal(), v) ≈ [1, 0] atol=atol rtol=rtol

        @test MOI.canget(m, MOI.ConstraintPrimal(), c)
        @test MOI.get(m, MOI.ConstraintPrimal(), c) ≈ 1 atol=atol rtol=rtol

        if MOI.get(solver, MOI.SupportsDuals())
            @test MOI.canget(m, MOI.DualStatus())
            @test MOI.get(m, MOI.DualStatus()) == MOI.FeasiblePoint
            @test MOI.canget(m, MOI.ConstraintDual(), c)
            @test MOI.get(m, MOI.ConstraintDual(), c) ≈ -1 atol=atol rtol=rtol

            # reduced costs
            @test MOI.canget(m, MOI.ConstraintDual(), vc1)
            @test MOI.get(m, MOI.ConstraintDual(), vc1) ≈ 0 atol=atol rtol=rtol
            @test MOI.canget(m, MOI.ConstraintDual(), vc2)
            @test MOI.get(m, MOI.ConstraintDual(), vc2) ≈ 1 atol=atol rtol=rtol
        end

        # change objective to Max +x

        objf = MOI.ScalarAffineFunction(v, [1.0,0.0], 0.0)
        MOI.set!(m, MOI.ObjectiveFunction(), objf)
        MOI.set!(m, MOI.ObjectiveSense(), MOI.MaxSense)

        if MOI.canget(m, MOI.ObjectiveFunction())
            @test objf ≈ MOI.get(m, MOI.ObjectiveFunction())
        end

        @test MOI.get(m, MOI.ObjectiveSense()) == MOI.MaxSense

        MOI.optimize!(m)

        @test MOI.canget(m, MOI.TerminationStatus())
        @test MOI.get(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.canget(m, MOI.PrimalStatus())
        @test MOI.get(m, MOI.PrimalStatus()) == MOI.FeasiblePoint

        @test MOI.canget(m, MOI.ObjectiveValue())
        @test MOI.get(m, MOI.ObjectiveValue()) ≈ 1 atol=atol rtol=rtol

        @test MOI.canget(m, MOI.VariablePrimal(), v)
        @test MOI.get(m, MOI.VariablePrimal(), v) ≈ [1, 0] atol=atol rtol=rtol

        if MOI.get(solver, MOI.SupportsDuals())
            @test MOI.canget(m, MOI.DualStatus())
            @test MOI.get(m, MOI.DualStatus()) == MOI.FeasiblePoint
            @test MOI.canget(m, MOI.ConstraintDual(), c)
            @test MOI.get(m, MOI.ConstraintDual(), c) ≈ -1 atol=atol rtol=rtol

            @test MOI.canget(m, MOI.ConstraintDual(), vc1)
            @test MOI.get(m, MOI.ConstraintDual(), vc1) ≈ 0 atol=atol rtol=rtol
            @test MOI.canget(m, MOI.ConstraintDual(), vc2)
            @test MOI.get(m, MOI.ConstraintDual(), vc2) ≈ 1 atol=atol rtol=rtol
        end

        # add new variable to get :
        # max x + 2z
        # s.t. x + y + z <= 1
        # x,y,z >= 0

        z = MOI.addvariable!(m)
        push!(v, z)
        @test v[3] == z

        vc3 = MOI.addconstraint!(m, MOI.SingleVariable(v[3]), MOI.GreaterThan(0.0))
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.SingleVariable,MOI.GreaterThan{Float64}}()) == 3

        @test MOI.canmodifyconstraint(m, c, MOI.ScalarCoefficientChange{Float64}(z, 1.0))
        MOI.modifyconstraint!(m, c, MOI.ScalarCoefficientChange{Float64}(z, 1.0))

        @test MOI.canmodifyobjective(m, MOI.ScalarCoefficientChange{Float64}(z, 2.0))
        MOI.modifyobjective!(m, MOI.ScalarCoefficientChange{Float64}(z, 2.0))

        @test MOI.get(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}}()) == 1
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.SingleVariable,MOI.GreaterThan{Float64}}()) == 3

        MOI.optimize!(m)

        @test MOI.canget(m, MOI.TerminationStatus())
        @test MOI.get(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.canget(m, MOI.ResultCount())
        @test MOI.get(m, MOI.ResultCount()) >= 1

        @test MOI.canget(m, MOI.PrimalStatus(1))
        @test MOI.get(m, MOI.PrimalStatus(1)) == MOI.FeasiblePoint

        @test MOI.canget(m, MOI.ObjectiveValue())
        @test MOI.get(m, MOI.ObjectiveValue()) ≈ 2 atol=atol rtol=rtol

        @test MOI.canget(m, MOI.VariablePrimal(), v)
        @test MOI.get(m, MOI.VariablePrimal(), v) ≈ [0, 0, 1] atol=atol rtol=rtol

        @test MOI.canget(m, MOI.ConstraintPrimal(), c)
        @test MOI.get(m, MOI.ConstraintPrimal(), c) ≈ 1 atol=atol rtol=rtol

        if MOI.get(solver, MOI.SupportsDuals())
            @test MOI.canget(m, MOI.DualStatus())
            @test MOI.get(m, MOI.DualStatus()) == MOI.FeasiblePoint
            @test MOI.canget(m, MOI.ConstraintDual(), c)
            @test MOI.get(m, MOI.ConstraintDual(), c) ≈ -2 atol=atol rtol=rtol

            @test MOI.canget(m, MOI.ConstraintDual(), vc1)
            @test MOI.get(m, MOI.ConstraintDual(), vc1) ≈ 1 atol=atol rtol=rtol
            @test MOI.canget(m, MOI.ConstraintDual(), vc2)
            @test MOI.get(m, MOI.ConstraintDual(), vc2) ≈ 2 atol=atol rtol=rtol
            @test MOI.canget(m, MOI.ConstraintDual(), vc3)
            @test MOI.get(m, MOI.ConstraintDual(), vc3) ≈ 0 atol=atol rtol=rtol
        end

        # setting lb of x to -1 to get :
        # max x + 2z
        # s.t. x + y + z <= 1
        # x >= -1
        # y,z >= 0
        @test MOI.canmodifyconstraint(m, vc1, MOI.GreaterThan(-1.0))
        MOI.modifyconstraint!(m, vc1, MOI.GreaterThan(-1.0))

        MOI.optimize!(m)

        @test MOI.canget(m, MOI.TerminationStatus())
        @test MOI.get(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.canget(m, MOI.ResultCount())
        @test MOI.get(m, MOI.ResultCount()) >= 1

        @test MOI.canget(m, MOI.PrimalStatus())
        @test MOI.get(m, MOI.PrimalStatus()) == MOI.FeasiblePoint

        @test MOI.canget(m, MOI.ObjectiveValue())
        @test MOI.get(m, MOI.ObjectiveValue()) ≈ 3 atol=atol rtol=rtol

        # put lb of x back to 0 and fix z to zero to get :
        # max x + 2z
        # s.t. x + y + z <= 1
        # x, y >= 0, z = 0
        @test MOI.canmodifyconstraint(m, vc1, MOI.GreaterThan(0.0))
        MOI.modifyconstraint!(m, vc1, MOI.GreaterThan(0.0))

        @test MOI.candelete(m, vc3)
        MOI.delete!(m, vc3)

        vc3 = MOI.addconstraint!(m, MOI.SingleVariable(v[3]), MOI.EqualTo(0.0))
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.SingleVariable,MOI.GreaterThan{Float64}}()) == 2

        MOI.optimize!(m)

        @test MOI.canget(m, MOI.TerminationStatus())
        @test MOI.get(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.canget(m, MOI.ResultCount())
        @test MOI.get(m, MOI.ResultCount()) >= 1

        @test MOI.canget(m, MOI.PrimalStatus())
        @test MOI.get(m, MOI.PrimalStatus()) == MOI.FeasiblePoint

        @test MOI.canget(m, MOI.ObjectiveValue())
        @test MOI.get(m, MOI.ObjectiveValue()) ≈ 1 atol=atol rtol=rtol

        # modify affine linear constraint set to be == 2 to get :
        # max x + 2z
        # s.t. x + y + z == 2
        # x,y >= 0, z = 0
        @test MOI.candelete(m, c)
        MOI.delete!(m, c)
        cf = MOI.ScalarAffineFunction(v, [1.0,1.0,1.0], 0.0)
        c = MOI.addconstraint!(m, cf, MOI.EqualTo(2.0))
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}}()) == 0
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}()) == 1

        MOI.optimize!(m)

        @test MOI.canget(m, MOI.TerminationStatus())
        @test MOI.get(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.canget(m, MOI.ResultCount())
        @test MOI.get(m, MOI.ResultCount()) >= 1

        @test MOI.canget(m, MOI.PrimalStatus())
        @test MOI.get(m, MOI.PrimalStatus()) == MOI.FeasiblePoint

        @test MOI.canget(m, MOI.ObjectiveValue())
        @test MOI.get(m, MOI.ObjectiveValue()) ≈ 2 atol=atol rtol=rtol

        # modify objective function to x + 2y to get :
        # max x + 2y
        # s.t. x + y + z == 2
        # x,y >= 0, z = 0

        objf = MOI.ScalarAffineFunction(v, [1.0,2.0,0.0], 0.0)
        MOI.set!(m, MOI.ObjectiveFunction(), objf)
        MOI.set!(m, MOI.ObjectiveSense(), MOI.MaxSense)

        MOI.optimize!(m)

        @test MOI.canget(m, MOI.TerminationStatus())
        @test MOI.get(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.canget(m, MOI.ResultCount())
        @test MOI.get(m, MOI.ResultCount()) >= 1

        @test MOI.canget(m, MOI.PrimalStatus())
        @test MOI.get(m, MOI.PrimalStatus()) == MOI.FeasiblePoint

        @test MOI.canget(m, MOI.ObjectiveValue())
        @test MOI.get(m, MOI.ObjectiveValue()) ≈ 4 atol=atol rtol=rtol

        @test MOI.canget(m, MOI.VariablePrimal(), v)
        @test MOI.get(m, MOI.VariablePrimal(), v) ≈ [0, 2, 0] atol=atol rtol=rtol

        # add constraint x - y >= 0 to get :
        # max x+2y
        # s.t. x + y + z == 2
        # x - y >= 0
        # x,y >= 0, z = 0

        cf2 = MOI.ScalarAffineFunction(v, [1.0, -1.0, 0.0], 0.0)
        c2 = MOI.addconstraint!(m, cf2, MOI.GreaterThan(0.0))
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.EqualTo{Float64}}()) == 1
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.GreaterThan{Float64}}()) == 1
        @test MOI.get(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}}()) == 0

        MOI.optimize!(m)

        @test MOI.canget(m, MOI.TerminationStatus())
        @test MOI.get(m, MOI.TerminationStatus()) == MOI.Success

        @test MOI.canget(m, MOI.ResultCount())
        @test MOI.get(m, MOI.ResultCount()) >= 1

        @test MOI.canget(m, MOI.PrimalStatus())
        @test MOI.get(m, MOI.PrimalStatus()) == MOI.FeasiblePoint

        @test MOI.canget(m, MOI.ObjectiveValue())
        @test MOI.get(m, MOI.ObjectiveValue()) ≈ 3 atol=atol rtol=rtol

        @test MOI.canget(m, MOI.VariablePrimal(), v)
        @test MOI.get(m, MOI.VariablePrimal(), v) ≈ [1, 1, 0] atol=atol rtol=rtol

        @test MOI.canget(m, MOI.ConstraintPrimal(), c)
        @test MOI.get(m, MOI.ConstraintPrimal(), c) ≈ 2 atol=atol rtol=rtol

        if MOI.get(solver, MOI.SupportsDuals())
            @test MOI.canget(m, MOI.DualStatus(1))
            @test MOI.get(m, MOI.DualStatus(1)) == MOI.FeasiblePoint

            @test MOI.canget(m, MOI.ConstraintDual(), c)
            @test MOI.get(m, MOI.ConstraintDual(), c) ≈ -1.5 atol=atol rtol=rtol
            @test MOI.get(m, MOI.ConstraintDual(), c2) ≈ 0.5 atol=atol rtol=rtol

            @test MOI.canget(m, MOI.ConstraintDual(), vc1)
            @test MOI.get(m, MOI.ConstraintDual(), vc1) ≈ 0 atol=atol rtol=rtol
            @test MOI.canget(m, MOI.ConstraintDual(), vc2)
            @test MOI.get(m, MOI.ConstraintDual(), vc2) ≈ 0 atol=atol rtol=rtol
            @test MOI.canget(m, MOI.ConstraintDual(), vc3)
            @test MOI.get(m, MOI.ConstraintDual(), vc3) ≈ 1.5 atol=atol rtol=rtol
        end
    end
end


linear1test(solver)
