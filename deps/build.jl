using BinaryProvider # requires BinaryProvider 0.3.0 or later

dependencies = [
     "https://github.com/JuliaOpt/OsiBuilder/releases/download/v0.107.9-1/build_OsiBuilder.v0.107.9.jl",
     "https://github.com/JuliaOpt/CoinUtilsBuilder/releases/download/v2.10.14-1/build_CoinUtilsBuilder.v2.10.14.jl",
     "https://github.com/JuliaOpt/COINGLPKBuilder/releases/download/v1.10.5-1/build_COINGLPKBuilder.v1.10.5.jl",
     "https://github.com/JuliaOpt/COINMumpsBuilder/releases/download/v1.6.0-1/build_COINMumpsBuilder.v1.6.0.jl",
     "https://github.com/JuliaOpt/COINMetisBuilder/releases/download/v1.3.5-1/build_COINMetisBuilder.v1.3.5.jl",
     "https://github.com/JuliaOpt/COINLapackBuilder/releases/download/v1.5.6-1/build_COINLapackBuilder.v1.5.6.jl",
     "https://github.com/JuliaOpt/COINBLASBuilder/releases/download/v1.4.6-1/build_COINBLASBuilder.v1.4.6.jl",
     "https://github.com/JuliaOpt/ASLBuilder/releases/download/v3.1.0-1/build_ASLBuilder.v3.1.0.jl"
 ]

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))
products = [
    LibraryProduct(prefix, ["libOsiClp"], :libOsiClp),
    LibraryProduct(prefix, ["libClp"], :libClp),
    LibraryProduct(prefix, ["libClpSolver"], :libClpSolver),
]

# Download binaries from hosted location
bin_prefix = "https://github.com/JuliaOpt/ClpBuilder/releases/download/v1.16.11-1"

# Listing of files generated by BinaryBuilder:
download_info = Dict(
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.aarch64-linux-gnu-gcc4.tar.gz", "7e23df1051244cc8bd63c86080d4be4ec93395a5287307f64fd3c2dfc5ff17af"),
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.aarch64-linux-gnu-gcc7.tar.gz", "9c5e663f6396fd4153b79b7c62bee1b5c3c6aa59ce2bd622ec04651b1f93b688"),
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.aarch64-linux-gnu-gcc8.tar.gz", "25ed2c904dd7ce2e04f58249b03470df41a1e780e9e681b7402cf3b2af9ba1b7"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.arm-linux-gnueabihf-gcc4.tar.gz", "13eba6ea26ef7b53400e18b5d8c477805242a3ebf52fb9fd281095a361f2b7f2"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.arm-linux-gnueabihf-gcc7.tar.gz", "70e2526aeb152643aa7252219ae0fa05ee3234ca4f3d02d6b0decc10a0d88a68"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.arm-linux-gnueabihf-gcc8.tar.gz", "c77b7cd2e098f06a2978e8a4d59d2e198e9725ccee7e75076f94c2ed0ddd88b3"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-linux-gnu-gcc4.tar.gz", "3de166e3fb6cfc5024a70e36412054ba42a418d6c6615c4b8181ca18de4713ca"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-linux-gnu-gcc7.tar.gz", "094e71013095e03643918511aa7e82c1d202c1ecf356957c6630bc19ec6230ed"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-linux-gnu-gcc8.tar.gz", "ab3e7e816cb0e9132e79826b20bbb692f671a26767fbb9ce0be641426fd740e7"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-w64-mingw32-gcc4.tar.gz", "9ffcf63f8c338e1ac557c7ee5d4524b9743934a4daad8d0f26e2aaf0e45b760a"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-w64-mingw32-gcc7.tar.gz", "bd19e5620944a3921f22c754dd81eaf61997610f980b18e21b32a64b3028a65a"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-w64-mingw32-gcc8.tar.gz", "8712972df47c1e56c8f8f00b371fa3d311debd8e8c365c7cd98c543074463bce"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-apple-darwin14-gcc4.tar.gz", "bcd454713c42f3f72d4e5a0d1bb1d4f98875eb4deea1730f84bc7159656a57ff"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-apple-darwin14-gcc7.tar.gz", "e59e70642664390a945016f84778dcccbd844958101106a12c0175583ea8fb9c"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-apple-darwin14-gcc8.tar.gz", "5dfa34bd4b01cc6988731d9f5dca896ef9965663408a4d7db656099ed52be11b"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-linux-gnu-gcc4.tar.gz", "6ba99b6edcce3390cf81f6968c78d2d1145b9828c0e66877216a987ed35a48a7"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-linux-gnu-gcc7.tar.gz", "e6054117c46b7972034963b23a74bde9502792dced9ec00ded185e48cce32515"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-linux-gnu-gcc8.tar.gz", "a6b15a10b86ed929078fc3cb75faac9d9b9e665f7fdcc18dcd0406d20fa8de8c"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-w64-mingw32-gcc4.tar.gz", "39aa0a1e9d7818e93ee5d162892df37c20c60a6169c7d1d7a5de59acd551963d"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-w64-mingw32-gcc7.tar.gz", "7f5e458cb389d6d043146cfdfd1b3301531f98fcfbbf343aa4bc1239879ec2ff"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-w64-mingw32-gcc8.tar.gz", "faa67966cedd814bcad547d3b65da9cc866e58facde95c6a55698f094fd065e5"),
)

# Install unsatisfied or updated dependencies:
unsatisfied = any(!satisfied(p; verbose=verbose) for p in products)
dl_info = choose_download(download_info, platform_key_abi())
if dl_info === nothing && unsatisfied
    # If we don't have a compatible .tar.gz to download, complain.
    # Alternatively, you could attempt to install from a separate provider,
    # build from source or something even more ambitious here.
    error("Your platform (\"$(Sys.MACHINE)\", parsed as \"$(triplet(platform_key_abi()))\") is not supported by this package!")
end

# If we have a download, and we are unsatisfied (or the version we're
# trying to install is not itself installed) then load it up!
if unsatisfied || !isinstalled(dl_info...; prefix=prefix)
    for dependency in reverse(dependencies)          # We do not check for already installed dependencies
       download(dependency,basename(dependency))
       evalfile(basename(dependency))
    end   
    # Download and install binaries
    install(dl_info...; prefix=prefix, force=true, verbose=verbose)
end

# Write out a deps.jl file that will contain mappings for our products
write_deps_file(joinpath(@__DIR__, "deps.jl"), products, verbose=verbose)
