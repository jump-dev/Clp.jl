println("WARNING: This build script has not been extensively tested. Please report any issues.")
clpname = "Clp-1.14.8"
prefix = joinpath(Pkg.dir(),"Clp","deps","usr")
libdir = joinpath(JULIA_HOME,"..","lib")
if !isfile("$clpname.tgz")
    run(`wget http://www.coin-or.org/download/source/Clp/$clpname.tgz`)
    run(`tar xvzf $clpname.tgz`)
    cd("$clpname")
    run(`cat ../Clp-makefile.patch` | `patch -p1`)
    run(`cat ../Clp-interface.patch` | `patch -p0`)
    run(`./configure --prefix=$prefix --with-blas="-L$libdir -lopenblas" --with-lapack=`)
    run(`make install`)
end

