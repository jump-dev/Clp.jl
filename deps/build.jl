using BinaryProvider # requires BinaryProvider 0.3.0 or later

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))
products = [
    LibraryProduct(prefix, ["libClp"], :libClp),
    LibraryProduct(prefix, ["libClpSolver"], :libClpSolver),
]

# Download binaries from hosted location
bin_prefix = "https://github.com/JuliaOpt/ClpBuilder/releases/download/v1.16.11-1a-static"

download_info = Dict(
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.aarch64-linux-gnu-gcc4.tar.gz", "9ab682ea624549801bf95c1169cd1da055941d9449b53185589e2c818d003e4b"),
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.aarch64-linux-gnu-gcc7.tar.gz", "92bd67f6f1ab11749b57228cd24134b5e8ca3a4b9ab695291528240ce33e6d68"),
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.aarch64-linux-gnu-gcc8.tar.gz", "7de4e0d406560ca36d68693371dc5d1aaafc8ccc92922ab2942a3ee17ed5cebb"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.arm-linux-gnueabihf-gcc4.tar.gz", "84ecf9257fd7cf35e0e03f70ba91330a120c202ca668bf4ea59ae462e114479d"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.arm-linux-gnueabihf-gcc7.tar.gz", "669c8eb49f1c55ae40f25b3102794771c8b85df4a01282f4de6e03a2d551dc6c"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.arm-linux-gnueabihf-gcc8.tar.gz", "1183075c704d653ba09e3f56971972ff8f8183a7081cca2aa8e244fc98c81be5"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-linux-gnu-gcc4.tar.gz", "5bb05730d63755d3c9a79dc8a78d09511568e62ad84ed802bc6f8077120c6b75"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-linux-gnu-gcc7.tar.gz", "a193103da7eee0b9c8ff7f331297d0e71baea40ef2cc25183e426f1f6b0e9e48"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-linux-gnu-gcc8.tar.gz", "5eed843e7db9459d1b47476499b313bfe9a2a020cc2fcb6cca5a6964529e12f1"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc6)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-w64-mingw32-gcc6.tar.gz", "a5c4036856d399bf06d33e422775a2e031ca906885af092690311a29c907dcee"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-w64-mingw32-gcc7.tar.gz", "e008d435230fa77da0897237df5b594c8f1a8828edbf324a9e1dc54f23db77d4"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-w64-mingw32-gcc8.tar.gz", "67c010fd535e371b44c5ad626107b3c42b24ae00f6f011cb8340cf1a68c8ee9b"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-apple-darwin14-gcc4.tar.gz", "07814295e81f24dd33dc1cfea349a840e086ba420808557161f3c7ada7b3cb03"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-apple-darwin14-gcc7.tar.gz", "a5259539475be446f6ade1c73513ebc8eb9d1bf7b8cd7f53430879ad3f2a8e4a"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-apple-darwin14-gcc8.tar.gz", "2ac3cd8ebc619887f2816b5f7cfdf99fe2d190ffe38eb92252d23d5ebfeb5715"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-linux-gnu-gcc4.tar.gz", "e2895e33c9193e4f547aaf8b75fe602c7506e8f87014dd9a92752c2a71cd6a63"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-linux-gnu-gcc7.tar.gz", "4c4f9c37c790a20d6e7874c684e898a012b5d07d61e125aa3a84fdf80620bcfc"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-linux-gnu-gcc8.tar.gz", "40dc02584678b264dde7be3084245a0b0efa84db4f143959fd1eca1a69c718a8"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc6)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-w64-mingw32-gcc6.tar.gz", "2ee820b182e4f878f973d66e6b8bea0a0cb9b3f9258de930bf7c94cbc80de71c"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-w64-mingw32-gcc7.tar.gz", "b4fe6a96a52ff7db307aab4cc6993f231f77ded6b39aeffa4b4c05cf11a8b770"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-w64-mingw32-gcc8.tar.gz", "2711d239bcbd3b7da18520694c447923c229c451b7f75c93573518f5424a6f61"),
)

# To fix gcc4 bug in Windows
# https://sourceforge.net/p/mingw-w64/bugs/727/
this_platform = platform_key_abi()
if typeof(this_platform)==Windows && this_platform.compiler_abi.gcc_version == :gcc4
   this_platform = Windows(arch(this_platform), libc=libc(this_platform), compiler_abi=CompilerABI(:gcc6))
end

# no dynamic dependencies until Pkg3 support for binaries
dependencies = [
#    "https://github.com/JuliaOpt/OsiBuilder/releases/download/v0.107.9-1/build_OsiBuilder.v0.107.9.jl",
#     "https://github.com/JuliaOpt/CoinUtilsBuilder/releases/download/v2.10.14-1/build_CoinUtilsBuilder.v2.10.14.jl",
#     "https://github.com/JuliaOpt/COINMumpsBuilder/releases/download/v1.6.0-1/build_COINMumpsBuilder.v1.6.0.jl",
#     "https://github.com/JuliaOpt/COINMetisBuilder/releases/download/v1.3.5-1/build_COINMetisBuilder.v1.3.5.jl",
#     "https://github.com/JuliaOpt/COINLapackBuilder/releases/download/v1.5.6-1/build_COINLapackBuilder.v1.5.6.jl",
#     "https://github.com/JuliaOpt/COINBLASBuilder/releases/download/v1.4.6-1/build_COINBLASBuilder.v1.4.6.jl",
#     "https://github.com/JuliaOpt/ASLBuilder/releases/download/v3.1.0-1/build_ASLBuilder.v3.1.0.jl"
]

custom_library = false
if haskey(ENV,"JULIA_CLP_LIBRARY_PATH")
    custom_products = [LibraryProduct(ENV["JULIA_CLP_LIBRARY_PATH"],product.libnames,product.variable_name) for product in products]
    if all(satisfied(p; verbose=verbose) for p in custom_products)
        products = custom_products
        custom_library = true
    else
        error("Could not install custom libraries from $(ENV["JULIA_CLP_LIBRARY_PATH"]).\nTo fall back to BinaryProvider call delete!(ENV,\"JULIA_CLP_LIBRARY_PATH\") and run build again.")
    end
end

if !custom_library
    # Install unsatisfied or updated dependencies:
    unsatisfied = any(!satisfied(p; verbose=verbose) for p in products)

    dl_info = choose_download(download_info, this_platform)
    if dl_info === nothing && unsatisfied
        # If we don't have a compatible .tar.gz to download, complain.
        # Alternatively, you could attempt to install from a separate provider,
        # build from source or something even more ambitious here.
        error("Your platform (\"$(Sys.MACHINE)\", parsed as \"$(triplet(platform_key_abi()))\") is not supported by this package!")
    end

    # If we have a download, and we are unsatisfied (or the version we're
    # trying to install is not itself installed) then load it up!
    if unsatisfied || !isinstalled(dl_info...; prefix=prefix)
        # Download and install binaries
        # no dynamic dependencies until Pkg3 support for binaries
        # for dependency in reverse(dependencies)          # We do not check for already installed dependencies
        #    download(dependency,basename(dependency))
        #    evalfile(basename(dependency))
        # end
        install(dl_info...; prefix=prefix, force=true, verbose=verbose)
    end
 end

# Write out a deps.jl file that will contain mappings for our products
write_deps_file(joinpath(@__DIR__, "deps.jl"), products, verbose=verbose)# using BinaryProvider # requires BinaryProvider 0.3.0 or later
