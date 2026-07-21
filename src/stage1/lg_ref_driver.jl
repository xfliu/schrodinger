# Double-precision LG reference at N=32,48,64 (N'=N), using certified rho = Stage-A inf(L2).
# Provides reference L1_LG values to validate the interval version's midpoint.
include("h2plus_lg_bounds.jl")
using Printf, JSON
using .H2plusLG
# certified inf(L2) from Stage-A certified_precision.json:
rho_cert = Dict(32=>-0.47230040183952, 48=>-0.40588252362138, 64=>-0.37792973588213)
configs = [(32,32,1.0),(48,48,1.0),(64,64,1.0)]
out=[]
for (N,Np,c0) in configs
    rec = H2plusLG.lg_bound(N,Np,c0; rho_override=rho_cert[N])
    @printf("N=%d Np=%d | L1_LG=%.10f rho=%.6f mu1N=%.10f B=%.4e nu=%.5f cg_res=%.2e (%.1fs)\n",
            N,Np,rec["L1_LG"],rec["rho"],rec["mu1N"],rec["B"],rec["nu"],rec["cg_residual"],rec["seconds"])
    push!(out, rec)
    open("lg_ref.json","w") do f; JSON.print(f,out); end
end
println("LG_REF_DONE")
