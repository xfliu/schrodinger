# dirichlet_lg.jl — certified interval Lehmann-Goerisch LOWER bounds for Dirichlet
# eigenvalues of H2+ on the truncated box Omega2 (sine basis, H^1_0).
#
# Sector structure (float, N=24): lambda1^D=-0.5469 is the ooo-sector ground;
# lambda2^D=-0.3275 is the eoo-sector ground (odd-in-x). Both are sector grounds =>
# each certifiable by single-test-vector LG within its own invariant subspace.
#
# Two-stage chain:
#  (A) certify_lambda2D  (eoo sector):
#        Galerkin lehmann_behnke r=1,2 -> mu_eoo1 (=lambda2^D), mu_eoo2 enclosures
#        Ceps_eoo (kinetic-scaled ground) + sharp_L -> separator for the eoo-ground LG
#        single-vector LG on eoo ground -> certified lambda2^D(Omega2) lower bound
#  (B) certify_lambda1D  (ooo sector), given rho = certified lambda2^D lower bound:
#        single-vector LG on ooo ground -> certified lambda1^D(Omega2) lower bound
#
# Goerisch A2 = <v,H'^{-1}v> via SPD identity (no interval inverse): CG solve against the
# SHIFTED mid-operator H'_mid = mid(H+cI); r = v - H'w; A2 = 2<w,v>-<w,H'w>+[0,||r||^2/lam_min].
module DirichletLG
using LinearAlgebra, Printf, SparseArrays
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
include("dirichlet_assembly.jl")
using .DirichletAssembly
const DA = DirichletAssembly
using .DirichletAssembly.MomentsVerified
using Arpack, IterativeSolvers
import Veigs: lehmann_behnke

const IV = Interval{Float64}
const LX = 10.0
# in-place: add (kin_scale*kin + sigma) to the diagonal of P (mutates P). Returns P as H.
function shift_diag!(P, kin; kin_scale=interval(1.0), sigma=interval(0.0))
    @inbounds for i in 1:length(kin)
        P[i,i] = P[i,i] + kin_scale*kin[i] + sigma
    end
    return P
end
nu_star(N::Int) = (interval(N+1)*interval(pi)/(interval(2.0)*interval(LX)))^2
fin(x)=isfinite(x) ? x : "nonfinite"

function sharp_L_iv(m1::IV, mk::IV, N::Int, eps::Float64, Ceps::IV)
    ε=interval(eps); one_=interval(1.0); σ=Ceps
    ms1=m1+σ; inf(ms1)>0 || return nothing
    etaV=sqrt(ε/(one_-ε)+Ceps/ms1); g=etaV/sqrt(ms1)
    Ch2=(one_+g*g)/((one_-ε)*nu_star(N)); msk=mk+σ
    return msk/(one_+Ch2*msk)-σ
end

function primary_block_indices(N::Int, Nprime::Int, px::Symbol, py::Symbol, pz::Symbol)
    m=N+1; mb=Nprime+1
    nxs=length(DA.sector_modes(m,px)); nys=length(DA.sector_modes(m,py)); nzs=length(DA.sector_modes(m,pz))
    nxb=length(DA.sector_modes(mb,px)); nyb=length(DA.sector_modes(mb,py))
    idx=Int[]
    for cz in 1:nzs, cy in 1:nys, cx in 1:nxs
        push!(idx,(cz-1)*nyb*nxb+(cy-1)*nxb+cx)
    end; idx
end

function lb_index(Hm, Hi, D, k, rho, sigma; nev=4)
    vals,vecs=eigs(Symmetric(Hm); nev=nev, which=:SR, maxiter=3000)
    p=sortperm(real(vals)); vals=real(vals[p]); vecs=real(vecs[:,p])
    V=reshape(vecs[:,k],D,1); Bi=interval.(sparse(1.0I,D,D))
    lg=lehmann_behnke(Hi,Bi,vals,V,rho,sigma,1.0,k,k; do_shift=true)
    return lg[1], vecs
end

# single-test-vector LG lower bound on the ground of sector (px,py,pz), separator rho_sep.
# Assembles at aux Nprime; test vector = primary-N sector ground, zero-padded.
# mu_lo_cert = rigorous lower bound on the finite matrix's smallest eigenvalue (for A2 lam_min).
function lg_sector_ground(Hc, Db, N, Nprime, px, py, pz; c, rho_sep, mu_lo_cert, v, Hms=nothing, verbose=true)
    # Hc is the interval H' = H + cI (already shifted in place); v is the (float) test vector.
    Hms === nothing && (Hms = mid.(Hc))
    vI=interval.(v)
    Hv=DA.imatvec(Hc,vI); A0=DA.idot(vI,Hv)
    w=zeros(Float64,Db); cg!(w, Symmetric(Hms), v; reltol=1e-12, maxiter=8000)
    cg_res=norm(Hms*w.-v)/norm(v)
    wI=interval.(w); Hw=DA.imatvec(Hc,wI); r=vI.-Hw
    wv=DA.idot(wI,vI); wHw=DA.idot(wI,Hw); r2=DA.idot(r,r)
    lam_min_lb=c+mu_lo_cert; @assert lam_min_lb>0 "shift too small: $lam_min_lb"
    corr=interval(0.0, sup(r2)/lam_min_lb); A2=interval(2.0)*wv-wHw+corr
    rho_p=interval(rho_sep)+interval(c)
    A=A0-rho_p; B=A0-interval(2.0)*rho_p+rho_p*rho_p*A2; nu=A/B
    L=(rho_p-rho_p/(interval(1.0)-nu))-interval(c)
    Bpos=inf(B)>0; nult1=sup(nu)<1
    verbose && @printf("  [LG %s%s%s] A0=[%.10f,%.10f] A2=[%.10f,%.10f] B=[%.3e,%.3e](+%s) nu<1=%s cg=%.1e => lower=%.12f w=%.2e\n",
        px,py,pz,inf(A0),sup(A0),inf(A2),sup(A2),inf(B),sup(B),Bpos,nult1,cg_res,inf(L),diam(L)); flush(stdout)
    return Dict{String,Any}("A0"=>[inf(A0),sup(A0)],"A2"=>[inf(A2),sup(A2)],
        "B"=>[fin(inf(B)),fin(sup(B))],"nu"=>[fin(inf(nu)),fin(sup(nu))],
        "B_positive"=>Bpos,"nu_lt_1"=>nult1,"cg_residual"=>cg_res,"rho_sep"=>rho_sep,
        "LG_lower"=>fin(inf(L)),"L_LG"=>[fin(inf(L)),fin(sup(L))],"L_LG_width"=>fin(diam(L)))
end

# Stage A: certify lambda2^D lower bound (eoo sector), returning the separator for stage B.
function certify_lambda2D(N::Int, Nprime::Int; c::Float64=1.0, rho2::Float64=-0.20, sigma2::Float64=-0.10,
        eps_grid=[0.40,0.30], nt::Int=48, npanel::Int=48,
        nodes::Vector{IV}, wts::Vector{IV}, verbose=true)
    T0=time()
    ta=@elapsed (P,kin,Db,dims)=DA.assemble_PK(Nprime,:even,:odd,:odd; nt=nt,npanel=npanel,nodes=nodes,wts=wts)
    verbose && @printf("  [eoo assemble N'=%d] D=%d dims=%s (%.1fs)\n",Nprime,Db,string(dims),ta); flush(stdout)
    # Ceps FIRST (needs kinetic-scaled matrices; do before mutating P into H).
    Ceps=Dict{Float64,IV}()
    for e in eps_grid
        Ha=DA.form_H(P,kin; kin_scale=interval(e)); Hma=mid.(Ha)
        eta,_=lb_index(Hma,Ha,Db,1,-1.6,-0.6); Ha=nothing; Hma=nothing; GC.gc()
        Ceps[e]=interval(-sup(eta),-inf(eta))
        verbose && @printf("  [eoo Ceps eps=%.2f] eta=[%.10f,%.10f] Ceps=[%.10f,%.10f]\n",
            e,inf(eta),sup(eta),inf(Ceps[e]),sup(Ceps[e])); flush(stdout)
    end
    # Build H in place on P (single interval matrix), Galerkin eoo ground (=lambda2^D) and eoo-2nd.
    shift_diag!(P, kin); Hm=mid.(P)                    # P is now H
    mu1N,_=lb_index(Hm,P,Db,1,-0.45,-0.28)
    mu2N,_=lb_index(Hm,P,Db,2,rho2,sigma2)
    verbose && @printf("  [eoo Galerkin] mu1=[%.12f,%.12f] mu2=[%.12f,%.12f]\n",
        inf(mu1N),sup(mu1N),inf(mu2N),sup(mu2N)); flush(stdout)
    # test vector from primary-N block of Hm (=mid H), before shifting to H'
    pidx=primary_block_indices(N,Nprime,:even,:odd,:odd)
    Hsub=Hm[pidx,pidx]; esub=eigen(Symmetric(Hsub)); vs=esub.vectors[:,1]
    v=zeros(Float64,Db); v[pidx].=vs; v./=norm(v)
    Hm=nothing; GC.gc()
    # sharp-L separator for the eoo-ground LG (base=mu1N, target=mu2N)
    bsep=nothing; be=NaN
    for e in eps_grid
        L=sharp_L_iv(mu1N,mu2N,Nprime,e,Ceps[e]); L===nothing && continue
        (bsep===nothing || inf(L)>inf(bsep)) && (bsep=L; be=e)
    end
    sep=inf(bsep); sep_valid = sep>sup(mu1N)
    verbose && @printf("  [eoo-2nd sharp-L] sep=%.12f (eps=%.2f) valid=%s\n",sep,be,sep_valid); flush(stdout)
    # shift P in place to H' = H + cI, then LG on eoo ground
    shift_diag!(P, kin; kin_scale=interval(0.0), sigma=interval(c))   # add c*I only (kin already added)
    lg=lg_sector_ground(P,Db,N,Nprime,:even,:odd,:odd; c=c,rho_sep=sep,mu_lo_cert=inf(mu1N),v=v,verbose=verbose)
    P=nothing; kin=nothing; GC.gc()
    merge!(lg,Dict{String,Any}("stage"=>"lambda2D","N"=>N,"Nprime"=>Nprime,"D_aux"=>Db,
        "mu1N_galerkin"=>[inf(mu1N),sup(mu1N)],"mu2N_galerkin"=>[inf(mu2N),sup(mu2N)],
        "eoo2nd_sharpL_sep"=>sep,"sep_eps"=>be,"sep_valid"=>sep_valid,
        "Ceps"=>Dict(string(e)=>[inf(Ceps[e]),sup(Ceps[e])] for e in eps_grid),
        "t_assemble"=>ta,"wall_total"=>time()-T0))
    return lg
end

# Stage B: certify lambda1^D lower bound (ooo sector) with given separator.
function certify_lambda1D(N::Int, Nprime::Int, rho_sep::Float64; c::Float64=1.0,
        mu_lo_floor::Float64=-0.65, nt::Int=48, npanel::Int=96, nodes::Vector{IV}, wts::Vector{IV}, verbose=true)
    T0=time()
    ta=@elapsed (P,kin,Db,dims)=DA.assemble_PK(Nprime,:odd,:odd,:odd; nt=nt,npanel=npanel,nodes=nodes,wts=wts)
    verbose && @printf("  [ooo assemble N'=%d] D=%d dims=%s (%.1fs)\n",Nprime,Db,string(dims),ta); flush(stdout)
    # single interval matrix: build H in place, primary-block Ritz vector, then shift in place to H'.
    shift_diag!(P, kin)                                # P is now H  (kinetic + Coulomb + 0 shift)
    pidx=primary_block_indices(N,Nprime,:odd,:odd,:odd)
    HmH=mid.(P)                                        # mid of H (for primary eigen only)
    Hsub=HmH[pidx,pidx]; esub=eigen(Symmetric(Hsub)); vs=esub.vectors[:,1]
    ritz_primary=esub.values[1]
    v=zeros(Float64,Db); v[pidx].=vs; v./=norm(v)
    HmH=nothing; GC.gc()
    verbose && @printf("  [ooo primary Ritz (float, block)] = %.12f\n",ritz_primary); flush(stdout)
    shift_diag!(P, kin; kin_scale=interval(0.0), sigma=interval(c))   # P is now H' = H + cI
    Hms=mid.(P)                                        # mid of H' (single copy, for CG)
    # rigorous mu_lo floor: Dirichlet lambda1^D >= Neumann lambda1 >= -0.6158 (Stage-A certified) => -0.65 safe
    lg=lg_sector_ground(P,Db,N,Nprime,:odd,:odd,:odd; c=c,rho_sep=rho_sep,mu_lo_cert=mu_lo_floor,v=v,Hms=Hms,verbose=verbose)
    P=nothing; Hms=nothing; kin=nothing; GC.gc()
    merge!(lg,Dict{String,Any}("stage"=>"lambda1D","N"=>N,"Nprime"=>Nprime,"D_aux"=>Db,
        "ritz_primary_float"=>ritz_primary,"mu_lo_floor"=>mu_lo_floor,"t_assemble"=>ta,"wall_total"=>time()-T0))
    return lg
end
end # module
