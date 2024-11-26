# Copyright (c) 2013: Clp.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module TestMOIWrapper

using Test

import Clp
import MathOptInterface as MOI

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_SolverName()
    @test MOI.get(Clp.Optimizer(), MOI.SolverName()) == "Clp"
    return
end

function test_supports_default_copy_to()
    @test !MOI.supports_incremental_interface(Clp.Optimizer())
    return
end

function test_runtests()
    # This is what JuMP would construct
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(Clp.Optimizer; with_bridge_type = Float64),
    )
    @test model.optimizer.model.model_cache isa
          MOI.Utilities.UniversalFallback{Clp.OptimizerCache}
    # `Variable.ZerosBridge` makes dual needed by some tests fail.
    MOI.Bridges.remove_bridge(
        model.optimizer,
        MOI.Bridges.Variable.ZerosBridge{Float64},
    )
    MOI.set(model, MOI.Silent(), true)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(; exclude = Any[MOI.ObjectiveBound]),
        exclude = [
            # Unable to prove infeasibility
            "test_conic_NormInfinityCone_INFEASIBLE",
            "test_conic_NormOneCone_INFEASIBLE",
        ],
    )
    return
end

function test_Nonexistant_unbounded_ray()
    inner = Clp.Optimizer()
    model =
        MOI.Utilities.CachingOptimizer(MOI.default_cache(inner, Float64), inner)
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
    return
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
    return
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
    return
end

function test_options_after_empty!()
    model = Clp.Optimizer()
    @test MOI.get(model, MOI.Silent()) == false
    MOI.set(model, MOI.Silent(), true)
    @test MOI.get(model, MOI.Silent()) == true
    MOI.empty!(model)
    @test MOI.get(model, MOI.Silent()) == true
    return
end

function test_attribute_TimeLimitSec()
    model = Clp.Optimizer()
    @test MOI.supports(model, MOI.TimeLimitSec())
    value = MOI.get(model, MOI.TimeLimitSec())
    MOI.set(model, MOI.TimeLimitSec(), 0.0)
    @test MOI.get(model, MOI.TimeLimitSec()) == 0.0
    MOI.set(model, MOI.TimeLimitSec(), nothing)
    @test MOI.get(model, MOI.TimeLimitSec()) === nothing
    MOI.set(model, MOI.TimeLimitSec(), 1.0)
    @test MOI.get(model, MOI.TimeLimitSec()) == 1.0
    MOI.set(model, MOI.TimeLimitSec(), value)
    return
end

function test_stopped_on_iterations()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    MOI.Utilities.loadfromstring!(
        inner,
        """
        variables: x1, x2, x3
        maxobjective: x1 + x2
        x1 + x2 + x3 <= 1.0
        x1 >= 0.0
        x2 >= 0.0
        x3 >= 0.0
        """
    )
    model = Clp.Optimizer()
    MOI.set(model, MOI.RawOptimizerAttribute("MaximumIterations"), 1)
    MOI.set(model, MOI.RawOptimizerAttribute("PresolveType"), 1)
    MOI.optimize!(model, inner)
    @test MOI.get(model, MOI.RawStatusString()) ==
          "3 - stopped on iterations etc"
    @test MOI.get(model, MOI.ResultCount()) == 0
    @test MOI.get(model, MOI.PrimalStatus()) == MOI.UNKNOWN_RESULT_STATUS
    @test MOI.get(model, MOI.DualStatus()) == MOI.UNKNOWN_RESULT_STATUS
    return
end

function test_isProvenPrimalInfeasible()
    inner = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    MOI.Utilities.loadfromstring!(
        inner,
        """
        variables: x, y
        minobjective: -1.0 * x + -1.0 * y
        -1.0 * x >= 1.0
        1.0 * y >= 1.0
        x >= 0.0
        y >= 0.0
        """
    )
    model = Clp.Optimizer()
    MOI.optimize!(model, inner)
    @test MOI.get(model, MOI.RawStatusString()) == "1 - primal infeasible"
    @test MOI.get(model, MOI.ResultCount()) == 1
    @test Clp.Clp_primalFeasible(model) == 0
    @test Clp.Clp_dualFeasible(model) == 0
    @test Clp.Clp_isProvenPrimalInfeasible(model) == 1
    return
end

end  # module TestMOIWrapper

TestMOIWrapper.runtests()
