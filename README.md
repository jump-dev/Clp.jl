Clp
=================

Interface to the **[Clp]** linear programming solver. Provides a complete interface to the low-level C API, as well as an implementation of the solver-independent ``MathProgSolverInterface`` type. The **[Cbc]** julia package is used to provide the binary dependencies; see that package's README for supported platforms and installation instructions. For users interested in a simple high-level ``linprog`` function, see the **[MathProgBase]** package. 

[Clp]: https://projects.coin-or.org/Clp
[Cbc]: https://github.com/mlubin/Cbc.jl
[MathProgBase]: https://github.com/JuliaOpt/MathProgBase.jl

### Using with MathProgBase

Clp provides a solver object that can be passed to ``linprog`` in MathProgBase (and used to create instances of the solver-independent ``AbstractMathProgModel`` type):

```julia
    using Clp
	using MathProgBase
	linprog(..., ClpSolver(Option1=value1,Option2=value2,...))
```

see the MathProgBase documentation for further information.

Options are solver-dependent. The following options are the most useful:

* ``PrimalTolerance`` - primal feasibility tolerance (default 1e-7)
* ``DualTolerance`` - dual feasibility tolerance (default 1e-7)
* ``DualObjectiveLimit`` - when using dual simplex (where the objective is monotonically changing), terminate when the objective exceeds this limit
* ``MaximumIterations`` - terminate after performing this number of simplex iterations
* ``MaximumSeconds`` - terminate after this many seconds have passed
* ``LogLevel`` - set to 1, 2, 3, or 4 for increasing output (default 0)
* ``PresolveType`` - set to 1 to disable presolve
* ``SolveType`` - choose the solution method

    - 0 - dual simplex
	- 1 - primal simplex
	- 3 - barrier with crossover to optimal basis
	- 4 - barrier without crossover to optimal basis
	- 5 - automatic

* ``InfeasibleReturn`` - set to 1 to return as soon as the problem is found to be infeasible (by default, an infeasibility proof is computed as well)

### Using the C interface

The low-level C interface is available in the ``ClpCInterface`` submodule:
```
    using Clp.ClpCInterface
```

Using this interface is only recommended for advanced users. The Julia API is essentially a thin wrapper around the interface exported by ``Clp/src/Clp_C_Interface.h``, which is documented in-line. 

