# Omega2 LG: aux N'=64 eee assembled once, primaries 32,48 from shared matrix.
using Printf, JSON, LinearAlgebra
include("lg_verified.jl")
using .LGVerified
using .LGVerified: IV, imatvec, idot, primary_block_indices
using .LGVerified.AssemblyVerified
const AV = LGVerified.AssemblyVerified
using .LGVerified.MomentsVerified
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
using Arpack, IterativeSolvers
const RHO_CERT = -0.54469513470285
const MU1_LO   = -0.64863173322543
const CSHIFT   = 1.0
Nprime = 64
prims  = [32,48]
nodes,wts = MomentsVerified.gl_reference(24)
println("GL reference certified (n_p=24)."); flush(stdout)
@printf("Assembling Omega2 aux eee N'=%d (npanel=96) ...\n", Nprime); flush(stdout)
ta = @elapsed (P,Kd,Db) = AV.assemble_PK(Nprime,:even,:even,:even; npanel=96, nodes=nodes, wts=wts)
@printf("  aux assemble N'=%d D=%d (%.1fs)\n", Nprime, Db, ta); flush(stdout)
Hc = AV.form_H(P, Kd; sigma=interval(CSHIFT)); P=nothing; Kd=nothing; GC.gc()
Hm = mid.(Hc)
out=[]
for N in prims
    T0=time()
    @printf("\n===== O2 LG primary N=%d aux N'=%d =====\n", N, Nprime); flush(stdout)
    pidx,_,_ = primary_block_indices(N, Nprime)
    Hsub = Hm[pidx,pidx]
    es = eigen(Symmetric(Hsub)); vs = es.vectors[:,1]
    v = zeros(Float64,Db); v[pidx].=vs; v ./= norm(v); vI = interval.(v)
    Hv = imatvec(Hc,vI); A0 = idot(vI,Hv)
    w = zeros(Float64,Db); cg!(w, Symmetric(Hm), v; reltol=1e-12, maxiter=8000)
    cg_res = norm(Hm*w .- v)/norm(v); wI = interval.(w)
    Hw = imatvec(Hc,wI); r = vI .- Hw
    wv=idot(wI,vI); wHw=idot(wI,Hw); r2=idot(r,r)
    lam_min_lb = CSHIFT + MU1_LO
    corr = interval(0.0, sup(r2)/lam_min_lb)
    A2 = interval(2.0)*wv - wHw + corr
    rho_p = interval(RHO_CERT) + interval(CSHIFT)
    A = A0 - rho_p
    B = A0 - interval(2.0)*rho_p + rho_p*rho_p*A2
    nu = A/B
    L1_LG = (rho_p - rho_p/(interval(1.0)-nu)) - interval(CSHIFT)
    Bpos=inf(B)>0; nult1=sup(nu)<1
    @printf("  A0=[%.12f,%.12f] A2=[%.12f,%.12f]\n",inf(A0),sup(A0),inf(A2),sup(A2))
    @printf("  B=[%.4e,%.4e](pos=%s) nu=[%.5f,%.5f](<1=%s) cg_res=%.2e\n",inf(B),sup(B),Bpos,inf(nu),sup(nu),nult1,cg_res)
    @printf("  L1_LG=[%.12f,%.12f] width=%.3e\n",inf(L1_LG),sup(L1_LG),diam(L1_LG)); flush(stdout)
    push!(out,Dict{String,Any}("N"=>N,"Nprime"=>Nprime,"c"=>CSHIFT,"D_aux"=>Db,
      "rho_L2cert"=>RHO_CERT,"mu1_lo_Nprime"=>MU1_LO,"lam_min_lb"=>lam_min_lb,
      "A0"=>[inf(A0),sup(A0)],"A2"=>[inf(A2),sup(A2)],"A2_corr_hi"=>sup(corr),
      "B"=>[inf(B),sup(B)],"nu"=>[inf(nu),sup(nu)],"B_positive"=>Bpos,"nu_lt_1"=>nult1,
      "cg_residual"=>cg_res,"r2_hi"=>sup(r2),
      "L1_LG"=>[inf(L1_LG),sup(L1_LG)],"L1_LG_lo"=>inf(L1_LG),"L1_LG_width"=>diam(L1_LG),
      "t_assemble"=>ta,"wall_lg"=>time()-T0))
    open("lg_verified_O2.json","w") do f; JSON.print(f,out,2); end
end
println("LG_O2_DONE")
