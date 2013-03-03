Clp
=================

Interface to the **[Clp]** linear programming solver. Provides a basic ``linprog`` interface (see in-line comments) as well as a complete interface to the low-level C API. This should work on Linux, OS X, and Windows, but it has not been widely tested. Binaries are provided for Windows (Vista-7-8). Windows users must install the Visual Studio **[redistributable]** package. OS X users will need a proper build environment. Please report any issues. 

[Clp]: https://projects.coin-or.org/Clp
[redistributable]: http://www.microsoft.com/en-us/download/details.aspx?id=30679

## To do

- Unified low-level interface to linear programming (like **[Osi]**)
- Improved linprog interface, allow some solver parameters and return more solution info

Contributions addressing any of these points (or others) are gladly accepted.

[Osi]: https://projects.coin-or.org/Osi
