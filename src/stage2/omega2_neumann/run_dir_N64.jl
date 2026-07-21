using Printf, JSON
include("dirichlet_verified.jl")
using .DirichletVerified
using .DirichletVerified.MomentsVerified
nodes,wts = MomentsVerified.gl_reference(24)
println("GL reference certified (n_p=24)."); flush(stdout)
@printf("\n===== Dirichlet N=64 Omega2 =====\n"); flush(stdout)
rec = DirichletVerified.dirichlet_upper(64; nodes=nodes, wts=wts)
open("dirichlet_N64.json","w") do f; JSON.print(f,[rec],2); end
@printf("  [checkpoint] N=64 Dirichlet UPPER=%.10f (%.1fs)\n",rec["dirichlet_upper_rigorous"],rec["wall_total"]); flush(stdout)
println("DIRICHLET_N64_DONE")
