using Printf, JSON
include("dirichlet_lg.jl"); using .DirichletLG
using .DirichletLG.MomentsVerified
nodes,wts=MomentsVerified.gl_reference(24)
RHO = -0.3382401620378125   # certified lambda2^D(Omega2) lower (Stage A)
println("===== Dirichlet lambda1^D lower (ooo), primary N=64 / aux N'=80 ====="); flush(stdout)
rB64=DirichletLG.certify_lambda1D(64,80, RHO; nodes=nodes,wts=wts)
open("dir_lambda1D_64.json","w") do f; JSON.print(f,[rB64],2); end
@printf("N64 LAMBDA1D_LOWER=%s Bpos=%s nu<1=%s\n",string(rB64["LG_lower"]),rB64["B_positive"],rB64["nu_lt_1"])
println("DIR_N64_DONE")
