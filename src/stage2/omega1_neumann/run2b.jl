using Printf, JSON
include("pipeline2b.jl")
using .Pipeline2b
using .Pipeline2b.MomentsVerified
Ns = length(ARGS)>0 ? parse.(Int,ARGS) : [32]
nodes,wts = MomentsVerified.gl_reference(24)
println("GL reference certified (n_p=24)."); flush(stdout)
results=Dict[]
for N in Ns
    @printf("\n========== N=%d ==========\n",N); flush(stdout)
    rec = Pipeline2b.run_N(N; nodes=nodes, wts=wts)
    push!(results, rec)
    open("certified_precision.json","w") do f; JSON.print(f,results,2); end
    @printf("  [checkpoint] N=%d (%.1fs)\n",N,rec["wall_total"]); flush(stdout)
end
println("ALL_N_DONE")
