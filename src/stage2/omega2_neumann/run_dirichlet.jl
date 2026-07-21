using Printf, JSON
include("dirichlet_verified.jl")
using .DirichletVerified
using .DirichletVerified.MomentsVerified
nodes,wts = MomentsVerified.gl_reference(24)
println("GL reference certified (n_p=24)."); flush(stdout)
Ns = length(ARGS)>0 ? parse.(Int,split(ARGS[1],",")) : [32]
out=[]
for N in Ns
    @printf("\n===== Dirichlet N=%d =====\n",N); flush(stdout)
    rec = DirichletVerified.dirichlet_upper(N; nodes=nodes, wts=wts)
    push!(out,rec)
    open("dirichlet_results.json","w") do f; JSON.print(f,out,2); end
    @printf("  [checkpoint] N=%d Dirichlet UPPER=%.10f (%.1fs)\n",N,rec["dirichlet_upper_rigorous"],rec["wall_total"]); flush(stdout)
end
println("DIRICHLET_DONE")
