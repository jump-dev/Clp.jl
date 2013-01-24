println("WARNING: This build script has not been extensively tested. Please report any issues.")
clpname = "Clp-1.14.8"
prefix = joinpath(Pkg.dir(),"Clp","deps","usr")
if !isfile("$clpname.tgz")
    run(`wget http://www.coin-or.org/download/source/Clp/$clpname.tgz`)
    run(`tar xvzf $clpname.tgz`)
    cd("$clpname")
    run(`cat ../Clp-makefile.patch` | `patch -p1`)
    run(`cat ../Clp-interface.patch` | `patch -p0`)
    run(`./configure --prefix=$prefix`)
    run(`make install`)
end

