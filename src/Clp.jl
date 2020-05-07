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

const _CLP_VERSION = VersionNumber(
    "$(Clp_VersionMajor()).$(Clp_VersionMinor()).$(Clp_VersionRelease())"
)

if v"1.17.6" <= _CLP_VERSION <= v"1.17.6"
    include("moi/MOI_wrapper.jl")
else
    @warn(
        "MOI wrapper not loaded because the version of Clp you have, " *
        "$_CLP_VERSION, is not supported by Clp.jl."
    )
end

end
