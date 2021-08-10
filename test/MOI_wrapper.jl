module TestMOIWrapper

using Test
using MathOptInterface
import Clp

const MOI = MathOptInterface

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

function test_SolverName()
    @test MOI.get(Clp.Optimizer(), MOI.SolverName()) == "Clp"
    return
end

function test_supports_default_copy_to()
    @test !MOI.supports_incremental_interface(Clp.Optimizer(), false)
    @test !MOI.supports_incremental_interface(Clp.Optimizer(), true)
    return
end

function test_runtests()
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            Clp.Optimizer(),
        ),
        Float64,
    )
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            exclude = Any[MOI.DualObjectiveValue, MOI.ObjectiveBound],
        ),
        exclude = [
            # TODO(odow): bug in Clp.jl
            "test_model_copy_to_UnsupportedAttribute",
            # Unable to prove infeasibility
            "test_conic_NormInfinityCone_INFEASIBLE",
            "test_conic_NormOneCone_INFEASIBLE",
        ],
    )
    return
end

function test_Nonexistant_unbounded_ray()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        Clp.Optimizer(),
    )
    MOI.set(model, MOI.Silent(), true)
    x = MOI.add_variables(model, 5)
    MOI.set(
        model,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, x), 0.0),
    )
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.optimize!(model)
    status = MOI.get(model, MOI.TerminationStatus())
    @test status == MOI.DUAL_INFEASIBLE
    return
end

function test_RawOptimizerAttribute()
    model = Clp.Optimizer()
    MOI.set(model, MOI.RawOptimizerAttribute("LogLevel"), 1)
    @test MOI.get(model, MOI.RawOptimizerAttribute("LogLevel")) == 1
    MOI.set(model, MOI.RawOptimizerAttribute("LogLevel"), 2)
    @test MOI.get(model, MOI.RawOptimizerAttribute("LogLevel")) == 2

    MOI.set(model, MOI.RawOptimizerAttribute("SolveType"), 1)
    @test MOI.get(model, MOI.RawOptimizerAttribute("SolveType")) == 1
    MOI.set(model, MOI.RawOptimizerAttribute("SolveType"), 4)
    @test MOI.get(model, MOI.RawOptimizerAttribute("SolveType")) == 4

    MOI.set(model, MOI.RawOptimizerAttribute("PresolveType"), 1)
    @test MOI.get(model, MOI.RawOptimizerAttribute("PresolveType")) == 1
    MOI.set(model, MOI.RawOptimizerAttribute("PresolveType"), 0)
    @test MOI.get(model, MOI.RawOptimizerAttribute("PresolveType")) == 0
end

function test_All_parameters()
    model = Clp.Optimizer()
    param = MOI.RawOptimizerAttribute("NotAnOption")
    @test !MOI.supports(model, param)
    @test_throws MOI.UnsupportedAttribute(param) MOI.get(model, param)
    @test_throws MOI.UnsupportedAttribute(param) MOI.set(model, param, false)
    for key in Clp.SUPPORTED_PARAMETERS
        @test MOI.supports(model, MOI.RawOptimizerAttribute(key))
        value = MOI.get(model, MOI.RawOptimizerAttribute(key))
        MOI.set(model, MOI.RawOptimizerAttribute(key), value)
        @test MOI.get(model, MOI.RawOptimizerAttribute(key)) == value
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

end  # module TestMOIWrapper

TestMOIWrapper.runtests()
