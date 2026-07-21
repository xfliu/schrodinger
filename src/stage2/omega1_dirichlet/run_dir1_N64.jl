using Printf, JSON
include("dirichlet_lg.jl"); using .DirichletLG
using .DirichletLG.MomentsVerified
nodes,wts=MomentsVerified.gl_reference(24)
RHO = -0.33122606818140476   # certified lambda2^D(Omega1) lower (Stage A)
println("===== Omega1 Stage B: lambda1^D lower (ooo), N=64/N'=80 ====="); flush(stdout)
rB64=DirichletLG.certify_lambda1D(64,80, RHO; nodes=nodes,wts=wts)
open("dir1_lambda1D_64.json","w") do f; JSON.print(f,[rB64],2); end
@printf("O1_N64_LAMBDA1D_LOWER=%s Bpos=%s nu<1=%s\n",string(rB64["LG_lower"]),rB64["B_positive"],rB64["nu_lt_1"])
println("O1_N64_DONE")
