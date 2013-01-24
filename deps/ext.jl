#@unix_only ENV["LD_LIBRARY_PATH"] = strcat(joinpath(Pkg.dir(),"Clp","deps","usr","lib"),":",ENV["LD_LIBRARY_PATH"])

const _jl_libClp = joinpath(Pkg.dir(),"Clp","deps","usr","lib","libClp.so")
