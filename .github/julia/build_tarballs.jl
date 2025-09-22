using BinaryBuilder, Pkg

sources = [
    GitSource(
        "https://github.com/coin-or/Clp.git",
        "7b9daa62d4c2710a368a17385913ce59d8c67b68",
    ),
]

script = raw"""
cd $WORKSPACE/srcdir/Clp*

# Remove misleading libtool files
rm -f ${prefix}/lib/*.la
update_configure_scripts

# old and custom autoconf
sed -i s/elf64ppc/elf64lppc/ configure

mkdir build
cd build/

export CPPFLAGS="${CPPFLAGS} -DNDEBUG -I${includedir} -I${includedir}/coin"
if [[ ${target} == *mingw* ]]; then
    export LDFLAGS="-L${libdir}"
elif [[ ${target} == *linux* ]]; then
    export LDFLAGS="-ldl -lrt"
fi

if [[ ${target} == *aarch64* ]] || [[ ${target} == *arm* ]]; then
    if [[ ${target} != *apple* ]]; then
        export CPPFLAGS="${CPPFLAGS} -D__arm__"
    fi
fi

# BLAS and LAPACK
if [[ "${target}" == *mingw* ]]; then
  LBT="-lblastrampoline-5"
else
  LBT="-lblastrampoline"
fi

../configure \
    --prefix=$prefix \
    --build=${MACHTYPE} \
    --host=${target} \
    --with-pic \
    --disable-pkg-config \
    --disable-debug \
    --disable-dependency-tracking \
    --enable-shared \
    lt_cv_deplibs_check_method=pass_all \
    --with-blas \
    --with-blas-lib="-L${libdir} ${LBT}" \
    --with-lapack \
    --with-lapack-lib="-L${libdir} ${LBT}" \
    --with-coinutils-lib="-lCoinUtils" \
    --with-osi-lib="-lCoinUtils -lOsi" \
    --with-mumps-lib="-L${libdir} -ldmumps -lmpiseq" \
    --with-mumps-incdir="${includedir}/libseq" \
    --with-metis-lib="-L${libdir} -lmetis" \
    --with-metis-incdir="${includedir}"

make -j${nproc}
make install

mv ${prefix}/share/coin/doc/* ${prefix}/share/licenses
"""

platforms = expand_cxxstring_abis(supported_platforms())
platforms = filter!(!Sys.isfreebsd, platforms)

products = [
    LibraryProduct("libClp", :libClp),
    LibraryProduct("libOsiClp", :libOsiClp),
    LibraryProduct("libClpSolver", :libClpSolver),
    ExecutableProduct("clp", :clp),
]

dependencies = [
    Dependency(PackageSpec(name="CoinUtils_jll", uuid="be027038-0da8-5614-b30d-e42594cb92df"), compat="=200.1100.1100"),
    Dependency(PackageSpec(name="METIS_jll", uuid="d00139f3-1899-568f-a2f0-47f597d42d70"), compat="=5.1.2"),
    Dependency(PackageSpec(name="MUMPS_seq_jll", uuid="d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d"), compat="=500.800.100"),
    Dependency(PackageSpec(name="Osi_jll", uuid="7da25872-d9ce-5375-a4d3-7a845f58efdd"), compat="=0.10800.1000"),
    Dependency(PackageSpec(name="libblastrampoline_jll", uuid="8e850b90-86db-534c-a0d3-1478176c7d93"), compat="5.4.0"),
    Dependency(PackageSpec(name="CompilerSupportLibraries_jll", uuid="e66e0078-7015-5450-92f7-15fbd957f2ae"))
]

build_tarballs(
    ARGS,
    "Clp",
    v"100.1700.901",
    sources,
    script,
    platforms,
    products,
    dependencies;
    preferred_gcc_version = v"8.1.0",
    preferred_llvm_version = v"13.0.1",
    julia_compat = "1.9",
)
