# COIN-OR Linear Programming Interface (Clp.jl)

![](https://www.coin-or.org/wordpress/wp-content/uploads/2014/08/COINOR.png)

`Clp.jl` is an interface to the **[COIN-OR Linear
Programming](https://projects.coin-or.org/Clp)** solver. It provides a complete
interface to the low-level C API, as well as an implementation of the
solver-independent `MathProgBase` and `MathOptInterface` API's.   

*Note: This wrapper is maintained by the JuliaOpt community and is not a COIN-OR
project.*

[![Build Status](https://travis-ci.org/JuliaOpt/Clp.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/Clp.jl)
[![codecov](https://codecov.io/gh/JuliaOpt/Clp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaOpt/Clp.jl)

[Clp]: https://projects.coin-or.org/Clp
[Cbc]: https://github.com/JuliaOpt/Cbc.jl

## Installation

The package is registered in `METADATA.jl` and so can be installed with `Pkg.add`.

```
julia> Pkg.add("Clp")
```

Clp.jl will use [BinaryProvider.jl](https://github.com/JuliaPackaging/BinaryProvider.jl) to automatically install the Clp binaries. This should work for both the official Julia binaries from `https://julialang.org/downloads/` and source-builds. 

## Custom Installation

To install custom built Clp binaries set the environmental variable `JULIA_CLP_LIBRARY_PATH` and call `import Pkg; Pkg.build("Clp")`. For instance, if the libraries are installed in `/opt/lib` just call
```julia
ENV["JULIA_Clp_LIBRARY_PATH"]="/opt/lib"
import Pkg; Pkg.build("Clp")
```
If you do not want BinaryProvider to download the default binaries on install set  `JULIA_CLP_LIBRARY_PATH`  before calling `import Pkg; Pkg.add("Clp")`.

To switch back to the default binaries clear `JULIA_CLP_LIBRARY_PATH` and call `import Pkg; Pkg.add("Clp")`.

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
