module Clp

import Clp_jll

function __init__()
    global libClp = Clp_jll.libClp
    version =
        VersionNumber("$(Clp_VersionMajor()).$(Clp_VersionMinor()).$(Clp_VersionRelease())")
    if !(v"1.17.2" <= version <= v"1.17.6")
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

for sym in names(@__MODULE__, all = true)
    sym_string = string(sym)
    if startswith(sym_string, "Clp_")
        @eval export $sym
    end
end

end
