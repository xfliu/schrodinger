# lg_oee.jl — certified lower bound on the TRUE second Neumann eigenvalue mu2(Omega)
# (oee-sector ground) via single-test-vector Lehmann-Goerisch on the oee block.
# The mu2 lower bound is the sharpest admissible separator rho for the eee-ground LG.
#
# One oee auxiliary assembly at N' supports:
#   lehmann_behnke r=s=1,2 -> verified Galerkin mu2^N', mu3^N' enclosures
#   sharp_L(mu2^N', mu3^N', Ceps) -> certified lower bound on TRUE oee-2nd
#         (= separator for the oee-ground LG)
#   single-vector LG on oee ground (primary-N block, zero-padded) -> mu2(Omega) lower bound.
module LGOEE
using LinearAlgebra, Printf, SparseArrays
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
include("assembly_verified.jl")
using .AssemblyVerified
const AV = AssemblyVerified
using .AssemblyVerified.MomentsVerified
using Arpack, IterativeSolvers
import Veigs: lehmann_behnke

const IV = Interval{Float64}
const LX = 20.0
nu_star(N::Int) = (interval(N+1)*interval(pi)/(interval(2.0)*interval(LX)))^2

function sharp_L_iv(m1::IV, mk::IV, N::Int, eps::Float64, Ceps::IV)
    ε=interval(eps); one_=interval(1.0); σ=Ceps
    ms1=m1+σ; inf(ms1)>0 || return nothing
    etaV=sqrt(ε/(one_-ε)+Ceps/ms1); g=etaV/sqrt(ms1)
    Ch2=(one_+g*g)/((one_-ε)*nu_star(N)); msk=mk+σ
    return msk/(one_+Ch2*msk)-σ
end

function primary_block_indices_oee(N::Int, Nprime::Int)
    nxs=length(AV.sector_idx(N,:odd)); nys=length(AV.sector_idx(N,:even)); nzs=nys
    nxb=length(AV.sector_idx(Nprime,:odd)); nyb=length(AV.sector_idx(Nprime,:even))
    idx=Int[]
    for cz in 1:nzs, cy in 1:nys, cx in 1:nxs
        push!(idx,(cz-1)*nyb*nxb+(cy-1)*nxb+cx)
    end; idx
end

function imatvec(M::Matrix{IV}, x::Vector{IV})
    D=length(x); y=Vector{IV}(undef,D)
    @inbounds for i in 1:D
        acc=interval(0.0); @simd for j in 1:D; acc+=M[i,j]*x[j]; end; y[i]=acc
    end; y
end
idot(a,b)=(s=interval(0.0);@inbounds for i in eachindex(a);s+=a[i]*b[i];end;s)

function lb_index(Hm, Hi, D, k, rho, sigma; nev=4)
    vals,vecs=eigs(Symmetric(Hm); nev=nev, which=:SR, maxiter=3000)
    p=sortperm(real(vals)); vals=real(vals[p]); vecs=real(vecs[:,p])
    V=reshape(vecs[:,k],D,1); Bi=interval.(sparse(1.0I,D,D))
    lg=lehmann_behnke(Hi,Bi,vals,V,rho,sigma,1.0,k,k; do_shift=true)
    return lg[1], vecs
end

function certify_mu2(N::Int, Nprime::Int; c::Float64=1.0,
        Ceps::Dict{Float64,IV}, mu2_lo_cert::Float64,
        eps_grid=[0.40,0.30], nt::Int=48, npanel::Int=96,
        nodes::Vector{IV}, wts::Vector{IV}, verbose=true)
    @assert Nprime>=N
    T0=time()
    ta=@elapsed (P,Kd,Db)=AV.assemble_PK_lean(Nprime,:odd,:even,:even; nt=nt,npanel=npanel,nodes=nodes,wts=wts)
    verbose && @printf("  [oee aux assemble N'=%d] D=%d (%.1fs)\n",Nprime,Db,ta); flush(stdout)
    Hoee=AV.form_H(P,Kd)
    Hm=mid.(Hoee)

    # verified Galerkin enclosures + Ritz vectors
    mu2N, vecs = lb_index(Hm, Hoee, Db, 1, -0.60, -0.25)
    mu3N, _    = lb_index(Hm, Hoee, Db, 2, -0.25, -0.10)
    verbose && @printf("  [oee Galerkin] mu2^N'=[%.12f,%.12f] mu3^N'=[%.12f,%.12f]\n",
        inf(mu2N),sup(mu2N),inf(mu3N),sup(mu3N)); flush(stdout)

    # certified lower bound on TRUE oee-2nd (separator for the oee-ground LG)
    bsep=nothing; be=NaN
    for e in eps_grid
        L=sharp_L_iv(mu2N,mu3N,Nprime,e,Ceps[e]); L===nothing && continue
        (bsep===nothing || inf(L)>inf(bsep)) && (bsep=L; be=e)
    end
    sep=inf(bsep)
    sep_valid = sep>sup(mu2N)
    verbose && @printf("  [oee-2nd sharp-L] lower=%.12f (eps=%.2f) valid_sep=%s\n",sep,be,sep_valid); flush(stdout)

    # test vector: primary-N oee block ground, zero-padded into aux
    pidx=primary_block_indices_oee(N,Nprime)
    Hsub=Hm[pidx,pidx]; es=eigen(Symmetric(Hsub)); vs=es.vectors[:,1]
    v=zeros(Float64,Db); v[pidx].=vs; v./=norm(v); vI=interval.(v)

    # shifted H' = H + cI (reuse P,Kd); CG MUST use the shifted mid-operator
    Hc=AV.form_H(P,Kd; sigma=interval(c)); P=nothing; Hoee=nothing; Hm=nothing; GC.gc()
    Hms=mid.(Hc)                                 # float mid of SHIFTED H' (SPD)

    Hv=imatvec(Hc,vI); A0=idot(vI,Hv)
    w=zeros(Float64,Db); cg!(w, Symmetric(Hms), v; reltol=1e-12, maxiter=8000)
    cg_res=norm(Hms*w.-v)/norm(v); Hms=nothing; GC.gc()
    wI=interval.(w); Hw=imatvec(Hc,wI); r=vI.-Hw
    wv=idot(wI,vI); wHw=idot(wI,Hw); r2=idot(r,r)
    lam_min_lb=c+mu2_lo_cert; @assert lam_min_lb>0
    corr=interval(0.0, sup(r2)/lam_min_lb); A2=interval(2.0)*wv-wHw+corr

    rho_p=interval(sep)+interval(c)
    A=A0-rho_p; B=A0-interval(2.0)*rho_p+rho_p*rho_p*A2; nu=A/B
    lam_hat=rho_p-rho_p/(interval(1.0)-nu); L2_LG=lam_hat-interval(c)
    Bpos=inf(B)>0; nult1=sup(nu)<1
    verbose && @printf("  [oee LG] A0=[%.12f,%.12f] A2=[%.12f,%.12f]\n",inf(A0),sup(A0),inf(A2),sup(A2))
    verbose && @printf("  [oee LG] B=[%.4e,%.4e](pos=%s) nu=[%.5f,%.5f](<1=%s) cg=%.2e\n",
        inf(B),sup(B),Bpos,inf(nu),sup(nu),nult1,cg_res)
    verbose && @printf("  [oee LG] mu2(Omega) >= %.12f  (L2_LG=[%.12f,%.12f] w=%.3e)\n",
        inf(L2_LG),inf(L2_LG),sup(L2_LG),diam(L2_LG)); flush(stdout)

    fin(x)=isfinite(x) ? x : "nonfinite"
    return Dict{String,Any}("N"=>N,"Nprime"=>Nprime,"c"=>c,"D_aux"=>Db,
        "mu2N_galerkin"=>[inf(mu2N),sup(mu2N)],"mu3N_galerkin"=>[inf(mu3N),sup(mu3N)],
        "oee2nd_sharpL_lower"=>sep,"sep_eps"=>be,"sep_valid"=>sep_valid,
        "A0"=>[inf(A0),sup(A0)],"A2"=>[inf(A2),sup(A2)],
        "B"=>[fin(inf(B)),fin(sup(B))],"nu"=>[fin(inf(nu)),fin(sup(nu))],"B_positive"=>Bpos,"nu_lt_1"=>nult1,
        "cg_residual"=>cg_res,"mu2_Omega_lower"=>fin(inf(L2_LG)),
        "L2_LG"=>[fin(inf(L2_LG)),fin(sup(L2_LG))],"L2_LG_width"=>fin(diam(L2_LG)),
        "t_assemble"=>ta,"wall_total"=>time()-T0)
end
end # module
