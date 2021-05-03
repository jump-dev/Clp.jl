function _fake_precompile()
    cache = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    clp = MOI.instantiate(Clp.Optimizer, with_bridge_type = Float64)
    MOI.copy_to(clp, cache; copy_names = false)
    return clp
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    T = Float64
    scalar_sets = (
        MOI.LessThan{T},
        MOI.GreaterThan{T},
        MOI.EqualTo{T},
        MOI.Interval{T},
    )
    scalar_functions = (MOI.SingleVariable, MOI.ScalarAffineFunction{T})
    MOI.precompile_model(
        MOI.Bridges.LazyBridgeOptimizer{MOI.Utilities.CachingOptimizer{
            Optimizer,
            MOI.Utilities.UniversalFallback{MOI.Utilities.Model{Float64}},
        }},
        [(F, S) for F in scalar_functions, S in scalar_sets],
    )
    @assert Base.precompile(_fake_precompile, ())
    return
end
