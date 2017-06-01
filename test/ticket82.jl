# See https://projects.coin-or.org/Clp/ticket/82
using Clp, Base.Test

m = MathProgBase.LinearQuadraticModel(ClpSolver())
MathProgBase.loadproblem!(m, "ticket82.mps")
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Unbounded
expunbray = [0,-1]
unbray = MathProgBase.getunboundedray(m)
@test !isempty(unbray) && norm(unbray - expunbray) < 1e-6
