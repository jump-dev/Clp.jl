# Copyright (c) 2013: Clp.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

function warmup()
    bridge = MOI.instantiate(Optimizer; with_bridge_type=Float64)
    no_bridge = MOI.instantiate(Optimizer)
    cache = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    MOI.optimize!(bridge, cache)
    MOI.optimize!(no_bridge, cache)
    return
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(warmup, ())
    return
end

_precompile_()
