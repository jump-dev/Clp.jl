# See https://projects.coin-or.org/Clp/ticket/79

m = MathProgBase.LinearQuadraticModel(ClpSolver())
MathProgBase.loadproblem!(m, "ticket79.mps")
MathProgBase.optimize!(m)
expinfray = [1,1,1,-1,-1,-1,0,0,0,0,0,0,0,0,0]
infray = MathProgBase.getinfeasibilityray(m)
@test !isempty(infray) && norm(infray - expinfray) < 1e-6
