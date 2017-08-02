using MathOptInterface, Clp, Base.Test
const MOI = MathOptInterface

solver = ClpMOISolver()

# @test MOI.supportsproblem(solver, MOI.ScalarAffineFunction{Float64}, [(MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}),(MOI.ScalarVariablewiseFunction,MOI.GreaterThan{Float64})])

m = MOI.SolverInstance(solver)

v = MOI.addvariables!(m, 2)
@test MOI.getattribute(m, MOI.NumberOfVariables()) == 2

cf = MOI.ScalarAffineFunction(v, [1.0,1.0], 0.0)
c = MOI.addconstraint!(m, cf, MOI.LessThan(1.0))

# @test MOI.getattribute(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.LessThan{Float64}}()) == 1

vc1 = MOI.addconstraint!(m, MOI.ScalarVariablewiseFunction(v[1]), MOI.GreaterThan(0.0))
vc2 = MOI.addconstraint!(m, MOI.ScalarVariablewiseFunction(v[2]), MOI.GreaterThan(0.0))

# @test MOI.getattribute(m, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64},MOI.GreaterThan{Float64}}()) == 2

objf = MOI.ScalarAffineFunction(v, [-1.0,0.0], 0.0)
MOI.setobjective!(m, MOI.MinSense, objf)

# @test MOI.getattribute(m, MOI.Sense()) == MOI.MinSense

@show Clp.ClpCInterface.get_constraint_matrix(m.inner)

# @show Clp.ClpCInterface.print_model(m.inner, "testlp")
