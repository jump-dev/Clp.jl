Clp
=================

Interface to the **[Clp]** linear programming solver. Provides a complete interface to the low-level C API, as well as an implementation of the solver-independent ``LinprogSolverInterface`` type. This package should work on Linux, OS X, and Windows, but it has not been widely tested. Binaries are provided for Windows (Vista-7-8). Windows users must install the Visual Studio **[redistributable]** package. OS X users will need a proper build environment. Please report any issues. The **[MathProgBase]** package provides a high-level ``linprog`` function.  

[Clp]: https://projects.coin-or.org/Clp
[redistributable]: http://www.microsoft.com/en-us/download/details.aspx?id=30679
[MathProgBase]: https://github.com/mlubin/MathProgBase.jl
