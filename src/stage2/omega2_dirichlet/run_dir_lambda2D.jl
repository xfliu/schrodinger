using Printf, JSON
include("dirichlet_lg.jl"); using .DirichletLG
using .DirichletLG.MomentsVerified
nodes,wts=MomentsVerified.gl_reference(24)
println("===== Dirichlet Stage A: certify lambda2^D(Omega2) lower bound (eoo, N=48/N'=64) ====="); flush(stdout)
rec=DirichletLG.certify_lambda2D(48,64; nodes=nodes,wts=wts)
open("dir_lambda2D_48_64.json","w") do f; JSON.print(f,[rec],2); end
@printf("LAMBDA2D_LOWER=%s sep_valid=%s Bpos=%s nu<1=%s\n",string(rec["LG_lower"]),rec["sep_valid"],rec["B_positive"],rec["nu_lt_1"])
println("DIR_STAGEA_DONE")
