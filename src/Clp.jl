# Copyright (c) 2013: Clp.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module Clp

import Clp_jll
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
    global libClp = Clp_jll.libClp
    version = VersionNumber(
        "$(Clp_VersionMajor()).$(Clp_VersionMinor()).$(Clp_VersionRelease())",
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
include("precompile.jl")

# Clp exports all `Clp_xxx` symbols. If you don't want all of these symbols in
# your environment, then use `import Clp` instead of `using Clp`.

for sym in filter(s -> startswith("$s", "Clp_"), names(@__MODULE__, all = true))
    @eval export $sym
end

end
