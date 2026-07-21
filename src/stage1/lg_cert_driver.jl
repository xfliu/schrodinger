
include("h2plus_lg_bounds.jl")
using Printf, JSON
using .H2plusLG
rho_cert = Dict(192=>-0.341979, 256=>-0.339859)
configs = [(192,256,1.0), (256,256,1.0)]
out = []
for (N,Np,c0) in configs
    rec_ritz = H2plusLG.lg_bound(N, Np, c0)
    rec_cert = H2plusLG.lg_bound(N, Np, c0; rho_override=rho_cert[N])
    @printf("N=%d Np=%d | rho=mu2:  L1_LG=%.7f (rho=%.6f B=%.3e nu=%.5f) | rho=L2cert: L1_LG=%.7f (rho=%.6f B=%.3e nu=%.5f)\n",
            N, Np, rec_ritz["L1_LG"], rec_ritz["rho"], get(rec_ritz,"B",NaN), get(rec_ritz,"nu",NaN),
            rec_cert["L1_LG"], rec_cert["rho"], get(rec_cert,"B",NaN), get(rec_cert,"nu",NaN))
    push!(out, Dict("N"=>N,"Nprime"=>Np,"ritz"=>rec_ritz,"cert"=>rec_cert))
end
open("lg_certified_rho.json","w") do f; JSON.print(f, out); end
println("WROTE lg_certified_rho.json")
