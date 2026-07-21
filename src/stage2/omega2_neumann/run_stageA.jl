using Printf, JSON
include("pipeline2b.jl")
using .Pipeline2b
using .Pipeline2b.MomentsVerified
nodes,wts = MomentsVerified.gl_reference(24)
println("GL reference certified (n_p=24)."); flush(stdout)
Ns = length(ARGS)>0 ? parse.(Int,split(ARGS[1],",")) : [32]
out=[]
for N in Ns
    @printf("\n===== Omega2 Stage-A N=%d =====\n",N); flush(stdout)
    rec = Pipeline2b.run_N(N; npanel=96, nodes=nodes, wts=wts)
    push!(out,rec)
    open("stageA_dom2.json","w") do f; JSON.print(f,out,2); end
    @printf("  [checkpoint] N=%d done bracket=[%.10f,%.10f]\n",N,rec["bracket_lo"],rec["bracket_hi"]); flush(stdout)
end
println("STAGEA_DOM2_DONE")
