using Printf, JSON
include("dirichlet_lg.jl"); using .DirichletLG
using .DirichletLG.MomentsVerified
nodes,wts=MomentsVerified.gl_reference(24)
println("===== Omega1 Stage A: certify lambda2^D(Omega1) lower (eoo, N=48/N'=64) ====="); flush(stdout)
rA=DirichletLG.certify_lambda2D(48,64; rho2=-0.15, sigma2=0.05, nodes=nodes,wts=wts)
open("dir1_lambda2D_48_64.json","w") do f; JSON.print(f,[rA],2); end
RHO = rA["LG_lower"]
@printf("O1_LAMBDA2D_LOWER=%s valid=%s Bpos=%s nu<1=%s\n",string(RHO),rA["sep_valid"],rA["B_positive"],rA["nu_lt_1"]); flush(stdout)
println("===== Omega1 Stage B: lambda1^D lower (ooo), N=48/N'=64 ====="); flush(stdout)
rB48=DirichletLG.certify_lambda1D(48,64, RHO; nodes=nodes,wts=wts)
open("dir1_lambda1D_48.json","w") do f; JSON.print(f,[rB48],2); end
@printf("O1_N48_LAMBDA1D_LOWER=%s Bpos=%s nu<1=%s\n",string(rB48["LG_lower"]),rB48["B_positive"],rB48["nu_lt_1"]); flush(stdout)
println("===== Omega1 Stage B: lambda1^D lower (ooo), N=64/N'=80 ====="); flush(stdout)
rB64=DirichletLG.certify_lambda1D(64,80, RHO; nodes=nodes,wts=wts)
open("dir1_lambda1D_64.json","w") do f; JSON.print(f,[rB64],2); end
@printf("O1_N64_LAMBDA1D_LOWER=%s Bpos=%s nu<1=%s\n",string(rB64["LG_lower"]),rB64["B_positive"],rB64["nu_lt_1"])
println("O1_ALL_DONE")
