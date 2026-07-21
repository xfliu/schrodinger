using Printf, JSON
include("lg_oee.jl"); using .LGOEE
using .LGOEE.MomentsVerified
using IntervalArithmetic: interval
nodes,wts=MomentsVerified.gl_reference(24)
println("GL reference certified (n_p=24)."); flush(stdout)
Ceps=Dict(0.40=>interval(0.8974058401931543,0.8974058402145533),
          0.30=>interval(1.0716666122207257,1.0716666122479594))
println("===== Omega2 oee LG: certify mu2 lower bound (N=64 / N'=64) ====="); flush(stdout)
rec=LGOEE.certify_mu2(64,64; Ceps=Ceps, mu2_lo_cert=-0.5446951347, nodes=nodes,wts=wts)
open("lg_oee_O2_N64.json","w") do f; JSON.print(f,[rec],2); end
@printf("MU2_OMEGA2_LOWER=%.12f sep_valid=%s Bpos=%s\n",rec["mu2_Omega_lower"],rec["sep_valid"],rec["B_positive"])
println("LG_OEE_DONE")
