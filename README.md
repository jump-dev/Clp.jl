Clp
=================

Interface to the **[Clp]** linear programming solver. Provides a complete interface to the low-level C API, as well as an implementation of the solver-independent ``MathProgSolverInterface`` type. The **[Cbc]** julia package is used to provide the binary dependencies; see that package's README for supported platforms and installation instructions. For users interested in a simple high-level ``linprog`` function, see the **[MathProgBase]** package.  

[Clp]: https://projects.coin-or.org/Clp
[Cbc]: https://github.com/mlubin/Cbc.jl
[MathProgBase]: https://github.com/mlubin/MathProgBase.jl
