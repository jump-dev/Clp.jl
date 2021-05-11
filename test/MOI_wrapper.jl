module TestMOIWrapper

using Test
using MathOptInterface
import Clp

const MOI = MathOptInterface

const OPTIMIZER = Clp.Optimizer()
MOI.set(OPTIMIZER, MOI.Silent(), true)

const CACHED = MOI.Utilities.CachingOptimizer(
    MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
    OPTIMIZER,
)

const BRIDGED = MOI.Bridges.full_bridge_optimizer(CACHED, Float64)

const CONFIG = MOI.Test.TestConfig(
    dual_objective_value = false,
    basis = true,
)

function test_SolverName()
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "Clp"
end

function test_supports_default_copy_to()
    @test !MOI.Utilities.supports_allocate_load(OPTIMIZER, false)
    @test !MOI.Utilities.supports_allocate_load(OPTIMIZER, true)
    @test !MOI.Utilities.supports_default_copy_to(OPTIMIZER, false)
    @test !MOI.Utilities.supports_default_copy_to(OPTIMIZER, true)
end

function test_basicconstraint()
    return MOI.Test.basic_constraint_tests(CACHED, CONFIG)
end

function test_unittest()
    return MOI.Test.unittest(
        BRIDGED,
        CONFIG,
        [
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
        ],
    )
end

function test_contlinear()
    return MOI.Test.contlineartest(BRIDGED, CONFIG, [
        # MOI.VariablePrimalStart not supported.
        "partial_start",
    ])
end

function test_nametest()
    return MOI.Test.nametest(BRIDGED)
end

function test_validtest()
    return MOI.Test.validtest(BRIDGED)
end

function test_emptytest()
    return MOI.Test.emptytest(BRIDGED)
end

function test_Nonexistant_unbounded_ray()
    MOI.empty!(BRIDGED)
    x = MOI.add_variables(BRIDGED, 5)
    MOI.set(
        BRIDGED,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0),
    )
    MOI.set(BRIDGED, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.optimize!(BRIDGED)
    status = MOI.get(BRIDGED, MOI.TerminationStatus())
    @test status == MOI.DUAL_INFEASIBLE
end

function test_RawParameter()
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

function test_All_parameters()
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

function test_copy_to_bug()
    model = MOI.Utilities.Model{Float64}()
    x = MOI.add_variable(model)
    con = [
        MOI.add_constraint(
            model,
            MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0),
            MOI.EqualTo(1.0),
        ) for i in 1:2
    ]
    clp = Clp.Optimizer()
    index_map = MOI.copy_to(clp, model)
    @test index_map[con[1]] != index_map[con[2]]
end

function test_options_after_empty!()
    model = Clp.Optimizer()
    @test MOI.get(model, MOI.Silent()) == false
    MOI.set(model, MOI.Silent(), true)
    @test MOI.get(model, MOI.Silent()) == true
    MOI.empty!(model)
    @test MOI.get(model, MOI.Silent()) == true
end

function test_farkas_dual_min()
    MOI.empty!(BRIDGED)
    model = BRIDGED
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.SingleVariable}(),
        MOI.SingleVariable(x[1]),
    )
    clb =
        MOI.add_constraint.(model, MOI.SingleVariable.(x), MOI.GreaterThan(0.0))
    c = MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([2.0, 1.0], x), 0.0),
        MOI.LessThan(-1.0),
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
    @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    clb_dual = MOI.get.(model, MOI.ConstraintDual(), clb)
    c_dual = MOI.get(model, MOI.ConstraintDual(), c)
    @show clb_dual, c_dual
    @test clb_dual[1] > 1e-6
    @test clb_dual[2] > 1e-6
    @test c_dual[1] < -1e-6
    @test clb_dual[1] ≈ -2 * c_dual atol = 1e-6
    @test clb_dual[2] ≈ -c_dual atol = 1e-6
end

function test_farkas_dual_min_interval()
    MOI.empty!(BRIDGED)
    model = BRIDGED
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.SingleVariable}(),
        MOI.SingleVariable(x[1]),
    )
    clb =
        MOI.add_constraint.(
            model,
            MOI.SingleVariable.(x),
            MOI.Interval(0.0, 10.0),
        )
    c = MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([2.0, 1.0], x), 0.0),
        MOI.LessThan(-1.0),
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
    @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    clb_dual = MOI.get.(model, MOI.ConstraintDual(), clb)
    c_dual = MOI.get(model, MOI.ConstraintDual(), c)
    @show clb_dual, c_dual
    @test clb_dual[1] > 1e-6
    @test clb_dual[2] > 1e-6
    @test c_dual[1] < -1e-6
    @test clb_dual[1] ≈ -2 * c_dual atol = 1e-6
    @test clb_dual[2] ≈ -c_dual atol = 1e-6
end

function test_farkas_dual_min_equalto()
    MOI.empty!(BRIDGED)
    model = BRIDGED
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.SingleVariable}(),
        MOI.SingleVariable(x[1]),
    )
    clb = MOI.add_constraint.(model, MOI.SingleVariable.(x), MOI.EqualTo(0.0))
    c = MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([2.0, 1.0], x), 0.0),
        MOI.LessThan(-1.0),
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
    @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    clb_dual = MOI.get.(model, MOI.ConstraintDual(), clb)
    c_dual = MOI.get(model, MOI.ConstraintDual(), c)
    @show clb_dual, c_dual
    @test clb_dual[1] > 1e-6
    @test clb_dual[2] > 1e-6
    @test c_dual[1] < -1e-6
    @test clb_dual[1] ≈ -2 * c_dual atol = 1e-6
    @test clb_dual[2] ≈ -c_dual atol = 1e-6
end

function test_farkas_dual_min_ii()
    MOI.empty!(BRIDGED)
    model = BRIDGED
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(-1.0, x[1])], 0.0),
    )
    clb = MOI.add_constraint.(model, MOI.SingleVariable.(x), MOI.LessThan(0.0))
    c = MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([-2.0, -1.0], x), 0.0),
        MOI.LessThan(-1.0),
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
    @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    clb_dual = MOI.get.(model, MOI.ConstraintDual(), clb)
    c_dual = MOI.get(model, MOI.ConstraintDual(), c)
    @show clb_dual, c_dual
    @test clb_dual[1] < -1e-6
    @test clb_dual[2] < -1e-6
    @test c_dual[1] < -1e-6
    @test clb_dual[1] ≈ 2 * c_dual atol = 1e-6
    @test clb_dual[2] ≈ c_dual atol = 1e-6
end

function test_farkas_dual_max()
    MOI.empty!(BRIDGED)
    model = BRIDGED
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.SingleVariable}(),
        MOI.SingleVariable(x[1]),
    )
    clb =
        MOI.add_constraint.(model, MOI.SingleVariable.(x), MOI.GreaterThan(0.0))
    c = MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([2.0, 1.0], x), 0.0),
        MOI.LessThan(-1.0),
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
    @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    clb_dual = MOI.get.(model, MOI.ConstraintDual(), clb)
    c_dual = MOI.get(model, MOI.ConstraintDual(), c)
    @show clb_dual, c_dual
    @test clb_dual[1] > 1e-6
    @test clb_dual[2] > 1e-6
    @test c_dual[1] < -1e-6
    @test clb_dual[1] ≈ -2 * c_dual atol = 1e-6
    @test clb_dual[2] ≈ -c_dual atol = 1e-6
end

function test_farkas_dual_max_ii()
    MOI.empty!(BRIDGED)
    model = BRIDGED
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 2)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(-1.0, x[1])], 0.0),
    )
    clb = MOI.add_constraint.(model, MOI.SingleVariable.(x), MOI.LessThan(0.0))
    c = MOI.add_constraint(
        model,
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([-2.0, -1.0], x), 0.0),
        MOI.LessThan(-1.0),
    )
    MOI.optimize!(model)
    @test MOI.get(model, MOI.TerminationStatus()) == MOI.INFEASIBLE
    @test MOI.get(model, MOI.DualStatus()) == MOI.INFEASIBILITY_CERTIFICATE
    clb_dual = MOI.get.(model, MOI.ConstraintDual(), clb)
    c_dual = MOI.get(model, MOI.ConstraintDual(), c)
    @show clb_dual, c_dual
    @test clb_dual[1] < -1e-6
    @test clb_dual[2] < -1e-6
    @test c_dual[1] < -1e-6
    @test clb_dual[1] ≈ 2 * c_dual atol = 1e-6
    @test clb_dual[2] ≈ c_dual atol = 1e-6
end

end  # module TestMOIWrapper

runtests(TestMOIWrapper)
