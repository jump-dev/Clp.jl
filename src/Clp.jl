# Copyright (c) 2013: Clp.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module Clp

import Clp_jll: libClp
import LinearAlgebra
import MathOptInterface as MOI
import OpenBLAS32_jll

function __init__()
    if VERSION >= v"1.9"
        config = LinearAlgebra.BLAS.lbt_get_config()
        if !any(lib -> lib.interface == :lp64, config.loaded_libs)
            LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)
        end
    end
    version = VersionNumber(
        Clp_VersionMajor(),
        Clp_VersionMinor(),
        Clp_VersionRelease(),
    )
    if !(v"1.17.2" <= version <= v"1.17.9")
        error(
            "You have installed version $version of Clp, which is not " *
            "supported by Clp.jl. If the version change was breaking, changes " *
            "will need to be made to the Julia code. Please open an issue at " *
            "https://github.com/jump-dev/Clp.jl.",
        )
    end
    return
end

include("libClp.jl")
include("MOI_wrapper/MOI_wrapper.jl")

# Clp exports all `Clp_xxx` symbols. If you don't want all of these symbols in
# your environment, then use `import Clp` instead of `using Clp`.

for sym in filter(s -> startswith("$s", "Clp_"), names(@__MODULE__, all = true))
    @eval export $sym
end

import PrecompileTools

PrecompileTools.@setup_workload begin
    PrecompileTools.@compile_workload begin
        let
            model = MOI.Utilities.CachingOptimizer(
                MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
                MOI.instantiate(Clp.Optimizer; with_bridge_type = Float64),
            )
            MOI.set(model, MOI.Silent(), true)
            x = MOI.add_variables(model, 3)
            sets = (MOI.GreaterThan(0.0), MOI.LessThan(2.0), MOI.EqualTo(1.0))
            for i in 1:3, f in (x[i], 1.0 * x[1] + 2.0 * x[2])
                MOI.supports_constraint(model, typeof(f), typeof(sets[i]))
                MOI.add_constraint(model, f, sets[i])
            end
            f = 1.0 * x[1] + x[2] + x[3]
            MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
            MOI.supports(model, MOI.ObjectiveFunction{typeof(f)}())
            MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
            MOI.optimize!(model)
            MOI.get(model, MOI.TerminationStatus())
            MOI.get(model, MOI.PrimalStatus())
            MOI.get(model, MOI.DualStatus())
            MOI.get(model, MOI.VariablePrimal(), x)
        end
    end
end

end
