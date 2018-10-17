# COIN-OR Linear Programming (Clp)

![](https://www.coin-or.org/wordpress/wp-content/uploads/2014/08/COINOR.png)

Interface to the **[COIN-OR Linear Programming](https://projects.coin-or.org/Clp)** solver. Provides a complete interface to the low-level C API, as well as an
implementation of the solver-independent `MathProgBase` and `MathOptInterface`
API's.   

[![Build Status](https://travis-ci.org/JuliaOpt/Clp.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/Clp.jl)

[![Clp](http://pkg.julialang.org/badges/Clp_0.6.svg)](http://pkg.julialang.org/?pkg=Clp&ver=0.6)

[Clp]: https://projects.coin-or.org/Clp
[Cbc]: https://github.com/JuliaOpt/Cbc.jl

## Installation

The package is registered in `METADATA.jl` and so can be installed with `Pkg.add`.

```
julia> Pkg.add("Clp")
```

Clp.jl will use [BinaryProvider.jl](https://github.com/JuliaPackaging/BinaryProvider.jl) to automatically install the Clp binaries.

## Custom Installation

After Clp.jl is installed and built, you can replace the installed binary dependencies with custom builds by overwritting the binaries and libraries in Clp.jl's `deps/usr` folder. For instance, Julia v0.6 this can be achieved by running
```bash
./configure --prefix=$HOME/.julia/v0.6/Clp/deps/usr ...
make
make install
```
in Clp's source folder.

Note that the custom binaries will not be overwritten by subsequent builds of the currently installed version of Clp.jl. However, if Clp.jl is updated and the update includes new BinaryProvider versions of the Clp binaries, then the custom binaries will be overwritten by the new BinaryProvider versions.

### Using with **[MathProgBase]**


Clp provides a solver object that can be passed to ``linprog`` in MathProgBase (and used to create instances of the solver-independent ``AbstractMathProgModel`` type):

```julia
using Clp
using MathProgBase
linprog(..., ClpSolver(Option1=value1,Option2=value2,...))
```

see the MathProgBase documentation for further information.

[MathProgBase]: https://github.com/JuliaOpt/MathProgBase.jl

Options are solver-dependent. The following options are the most useful (and well documented):

* ``PrimalTolerance`` - primal feasibility tolerance (default 1e-7)
* ``DualTolerance`` - dual feasibility tolerance (default 1e-7)
* ``DualObjectiveLimit`` - when using dual simplex (where the objective is monotonically changing), terminate when the objective exceeds this limit
* ``MaximumIterations`` - terminate after performing this number of simplex iterations
* ``MaximumSeconds`` - terminate after this many seconds have passed
* ``LogLevel`` - set to 1, 2, 3, or 4 for increasing output (default 0)
* ``PresolveType`` - set to 1 to disable presolve
* ``SolveType`` - choose the solution method:

    - 0 - dual simplex
	- 1 - primal simplex
	- 3 - barrier with crossover to optimal basis
	- 4 - barrier without crossover to optimal basis
	- 5 - automatic

* ``InfeasibleReturn`` - set to 1 to return as soon as the problem is found to be infeasible (by default, an infeasibility proof is computed as well)

### Using the C interface

The low-level C interface is available in the ``ClpCInterface`` submodule:
```julia
using Clp.ClpCInterface
```

Using this interface is only recommended for advanced users. The Julia API is essentially a thin wrapper around the interface exported by ``Clp/src/Clp_C_Interface.h``, which is documented in-line.
