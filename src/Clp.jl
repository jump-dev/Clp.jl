module Clp

if haskey(ENV,"JULIA_CLP_LIBRARY_PATH") || VERSION < v"1.3"
    if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
        include("../deps/deps.jl")
    else
        error("Clp not properly installed. Please run import Pkg; Pkg.build(\"Clp\")")
    end
else
    import Clp_jll: libClp
end

using CEnum

include("gen/ctypes.jl")
include("gen/libclp_common.jl")
include("gen/libclp_api.jl")

import MathOptInterface
const MOI = MathOptInterface
include("moi/MOI_wrapper.jl")

end
