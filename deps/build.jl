using BinDeps

@BinDeps.setup

@unix_only libclp = library_dependency("libclp",aliases=["libClp"])
@windows_only libclp = library_dependency("libclp",aliases=["CoinMP"])

@BinDeps.if_install begin

clpname = "Clp-1.14.8"

provides(Sources, URI("http://www.coin-or.org/download/source/Clp/$clpname.tgz"), 
    libclp, os = :Unix)

provides(Binaries, URI("http://www.mit.edu/~mlubin/CoinMP_julia.tar.gz"),
    libclp, os = :Windows)


prefix=joinpath(BinDeps.depsdir(libclp),"usr")
patchdir=BinDeps.depsdir(libclp)
clpsrcdir = joinpath(BinDeps.depsdir(libclp),"src",clpname) 

provides(SimpleBuild,
    (@build_steps begin
        #BinDeps.DirectoryRule(clpsrcdir,GetSources(libclp))
        GetSources(libclp)
        @build_steps begin
            ChangeDirectory(clpsrcdir)
            `cat $patchdir/Clp-makefile.patch` |> `patch -p1`
            `cat $patchdir/Clp-interface.patch` |> `patch -p0`
            `./configure --prefix=$prefix`
            `make install`
        end
    end),libclp, os = :Unix)

# TODO: Fix this
@windows_only begin
    # see CoinMP package for patches used to build this
    if !isfile("CoinMP_julia.tar.gz")
        run(download_cmd("http://www.mit.edu/~mlubin/CoinMP_julia.tar.gz","CoinMP_julia.tar.gz"))
    end
    if !isfile("CoinMP.dll")
        run(unpack_cmd("CoinMP_julia.tar.gz","."))
    end
end

@BinDeps.install

end

