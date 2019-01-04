using BinaryProvider # requires BinaryProvider 0.3.0 or later

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS
const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))
products = [
    LibraryProduct(prefix, ["libClp"], :libClp),
    LibraryProduct(prefix, ["libClpSolver"], :libClpSolver),
]

# Download binaries from hosted location
bin_prefix = "https://github.com/JuliaOpt/ClpBuilder/releases/download/v1.16.11-1-static"

download_info = Dict(
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.aarch64-linux-gnu-gcc4.tar.gz", "7486e8fee3ee02770276692784cbb2f6f5833445af5873ba4ee8e0dccce73ea9"),
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.aarch64-linux-gnu-gcc7.tar.gz", "d245f665a2cd764cb220a15bd25b8f22d39e29c1a8602830d76bccb98d6efcfa"),
    Linux(:aarch64, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.aarch64-linux-gnu-gcc8.tar.gz", "4f35d9c124ce56a7289c7648a417f1ecaa771282f66eb7f958ed5e4d434a9e0e"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.arm-linux-gnueabihf-gcc4.tar.gz", "476be15b10ee70b20f20ae4f5120e1664715de406bf504b8ac1ae447214aedaf"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.arm-linux-gnueabihf-gcc7.tar.gz", "ffdb45fd798892b850d95b2a26c59f3f86a702a82c0ab227a9f8aebb75b892f4"),
    Linux(:armv7l, libc=:glibc, call_abi=:eabihf, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.arm-linux-gnueabihf-gcc8.tar.gz", "fffbea4035857537006f3e466ef72369d4daa149e960bd69ee8c7d52562013f8"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-linux-gnu-gcc4.tar.gz", "cbcd28e3ac7a6bd1250294a804114181c87f1b616113225ef7d35dcfb7f6ffa7"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-linux-gnu-gcc7.tar.gz", "977b830cd6997272cdee726c6deab24a5ccb54e5d10626c47d6be630a684db44"),
    Linux(:i686, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-linux-gnu-gcc8.tar.gz", "718aad4ef7d46106ca8f1f082bface33c3b64e223b523ce02ba143e8e08b72e7"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc6)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-w64-mingw32-gcc6.tar.gz", "181fef44a57723a3d9e4cac6828aff308d2fa849245f6bf800110c7645a7bfdd"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-w64-mingw32-gcc7.tar.gz", "252ecd181bb32b5610d3131b3eee0c70cc68ae3f1d4fe2ba689256caf52afe86"),
    Windows(:i686, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.i686-w64-mingw32-gcc8.tar.gz", "0ac76950c0874f06bb52b13a85ec813181aa9c9162d741c8dc0d4bae01a22a26"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-apple-darwin14-gcc4.tar.gz", "5a9aeb6e27ff9677464bcf3ecc0735a67d5473937b8555b6dc6b9064ed921b82"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-apple-darwin14-gcc7.tar.gz", "262e7a7f16d08edf17d838a232016be29c6b5fc8484aa12790da2c44bed3917e"),
    MacOS(:x86_64, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-apple-darwin14-gcc8.tar.gz", "1ba927665b089b0e6738d5cc62eb8fa916e0d5387e38964b388205d44f10005b"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc4)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-linux-gnu-gcc4.tar.gz", "b7e1b16b7a478762cf24f3112f6d064d497b845bac63896fdda19085d2109985"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-linux-gnu-gcc7.tar.gz", "b955aadca6134586babe120c499aee6583715d4964075bf73bbc42216e5ba911"),
    Linux(:x86_64, libc=:glibc, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-linux-gnu-gcc8.tar.gz", "9b9ae725558a642b6a03a327236d91b1fc8c2c99c3afed932a3921d27ef1cbf6"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc6)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-w64-mingw32-gcc6.tar.gz", "03252a30570e9740dcc5f1413dd3481a68e8c12c7373fc1ac2a58cec73f3ff6f"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc7)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-w64-mingw32-gcc7.tar.gz", "76078cbbf17456ea3b27145a3539a1f16c26ac261a13485068c8718d55aba6a3"),
    Windows(:x86_64, compiler_abi=CompilerABI(:gcc8)) => ("$bin_prefix/ClpBuilder.v1.16.11.x86_64-w64-mingw32-gcc8.tar.gz", "705cd3fc170d1607899e9bcfe981c84ef5244a94187abbc479251c7f37b316e1"),
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
