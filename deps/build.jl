using BinDeps

@unix_only begin
clpname = "Clp-1.14.8"
prefix = joinpath(Pkg.dir(),"Clp","deps","usr")
libdir = joinpath(JULIA_HOME,"..","lib")
if !isfile("$clpname.tgz")
    run(download_cmd("http://www.coin-or.org/download/source/Clp/$clpname.tgz","$clpname.tgz"))
    run(`tar xvzf $clpname.tgz`)
    cd("$clpname")
    run(`cat ../Clp-makefile.patch` | `patch -p1`)
    run(`cat ../Clp-interface.patch` | `patch -p0`)
    # We should use Julia's blas/lapack, but this seems to cause some crashes
    #run(`./configure --prefix=$prefix --with-blas="-L$libdir -lopenblas" --with-lapack=`)
    run(`./configure --prefix=$prefix`)
    run(`make install`)
end
end # unix_only

@windows_only begin
    # see CoinMP package for patches used to build this
    if !isfile("CoinMP_julia.tar.gz")
        run(download_cmd("http://www.mit.edu/~mlubin/CoinMP_julia.tar.gz","CoinMP_julia.tar.gz"))
    end
    if !isfile("CoinMP.dll")
        run(unpack_cmd("CoinMP_julia.tar.gz","."))
    end
end


