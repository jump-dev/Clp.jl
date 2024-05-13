![](https://www.coin-or.org/wordpress/wp-content/uploads/2014/08/COINOR.png)

# Clp.jl

[![Build Status](https://github.com/jump-dev/Clp.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/Clp.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/Clp.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/Clp.jl)

[Clp.jl](https://github.com/jump-dev/Clp.jl) is a wrapper for the
[COIN-OR Linear Programming](https://projects.coin-or.org/Clp) solver.

The wrapper has two components:

 * a thin wrapper around the complete C API
 * an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

## Affiliation

This wrapper is maintained by the JuMP community and is not a COIN-OR project.

## Getting help

If you need help, please ask a question on the [JuMP community forum](https://jump.dev/forum).

If you have a reproducible example of a bug, please [open a GitHub issue](https://github.com/jump-dev/Clp.jl/issues/new).

## License

`Clp.jl` is licensed under the [MIT License](https://github.com/jump-dev/Clp.jl/blob/master/LICENSE.md).

The underlying solver, [coin-or/Clp](https://github.com/coin-or/Clp), is
licensed under the [Eclipse public license](https://github.com/coin-or/Clp/blob/master/LICENSE).

## Installation

Install Clp using `Pkg.add`:
```julia
import Pkg
Pkg.add("Clp")
```

In addition to installing the Clp.jl package, this will also download and
install the Clp binaries. You do not need to install Clp separately.

To use a custom binary, read the [Custom solver binaries](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
section of the JuMP documentation.

## Use with JuMP

To use Clp with JuMP, use `Clp.Optimizer`:
```julia
using JuMP, Clp
model = Model(Clp.Optimizer)
set_attribute(model, "LogLevel", 1)
set_attribute(model, "Algorithm", 4)
```

## MathOptInterface API

The Clp optimizer supports the following constraints and attributes.

List of supported objective functions:

 * [`MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}`](@ref)

List of supported variable types:

 * [`MOI.Reals`](@ref)

List of supported constraint types:

 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.Interval{Float64}`](@ref)
 * [`MOI.ScalarAffineFunction{Float64}`](@ref) in [`MOI.LessThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.EqualTo{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.GreaterThan{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.Interval{Float64}`](@ref)
 * [`MOI.VariableIndex`](@ref) in [`MOI.LessThan{Float64}`](@ref)

List of supported model attributes:

 * [`MOI.ObjectiveSense()`](@ref)

## Options

Options are, unfortunately, not well documented.

The following options are likely to be the most useful:

| Parameter            | Example | Explanation                                 |
| -------------------- | ------- | ------------------------------------------- |
| `PrimalTolerance`    | `1e-7`  | Primal feasibility tolerance                |
| `DualTolerance`      | `1e-7`  | Dual feasibility tolerance                  |
| `DualObjectiveLimit` | `1e308` | When using dual simplex (where the objective is monotonically changing), terminate when the objective exceeds this limit |
| `MaximumIterations`  | `2147483647` | Terminate after performing this number of simplex iterations |
| `MaximumSeconds`     | `-1.0`  | Terminate after this many seconds have passed. A negative value means no time limit |
| `LogLevel`           | `1`     | Set to 1, 2, 3, or 4 for increasing output. Set to `0` to disable output |
| `PresolveType`       | `0`     | Set to `1` to disable presolve              |
| `SolveType`          | `5`     | Solution method: dual simplex (`0`), primal simplex (`1`), sprint (`2`), barrier with crossover (`3`), barrier without crossover (`4`), automatic (`5`) |
| `InfeasibleReturn`   | `0`     | Set to 1 to return as soon as the problem is found to be infeasible (by default, an infeasibility proof is computed as well) |
| `Scaling`            | `3`     | `0` -off, `1` equilibrium, `2` geometric, `3` auto, `4` dynamic(later) |
| `Perturbation`       | `100`   | switch on perturbation (`50`), automatic (`100`), don't try perturbing (`102`) |

## C API

The C API can be accessed via `Clp.Clp_XXX` functions, where the names and
arguments are identical to the C API.
