# Double-precision LG reference with ENRICHED auxiliary space N'>N.
# This is where LG becomes non-degenerate and tightens below the Ritz value.
include("h2plus_lg_bounds.jl")
using Printf, JSON
using .H2plusLG
rho_cert = Dict(32=>-0.47230040183952, 48=>-0.40588252362138, 64=>-0.37792973588213)
# primary N, auxiliary N', shift c
configs = [(32,48,1.0),(32,64,1.0),(48,64,1.0)]
out=[]
for (N,Np,c0) in configs
    rec = H2plusLG.lg_bound(N,Np,c0; rho_override=rho_cert[N])
    @printf("N=%d Np=%d | L1_LG=%.10f (vs mu1N=%.10f, gap=%.3e) rho=%.6f B=%.4e nu=%.5f A2=%.6f cg_res=%.1e (%.1fs)\n",
            N,Np,rec["L1_LG"],rec["mu1N"],rec["mu1N"]-rec["L1_LG"],rec["rho"],rec["B"],rec["nu"],rec["A2"],rec["cg_residual"],rec["seconds"])
    push!(out, rec)
    open("lg_ref2.json","w") do f; JSON.print(f,out); end
end
println("LG_REF2_DONE")
