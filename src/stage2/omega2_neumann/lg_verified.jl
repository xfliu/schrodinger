# lg_verified.jl — certified interval-arithmetic single-test-vector Lehmann–Goerisch (LG)
# lower bound for the H2+ ground energy lambda_1. Rigorous counterpart of Stage-1
# h2plus_lg_bounds.jl. Reuses the validated interval assembly (assembly_verified.jl).
#
# Method (single normalized test vector v1 => LG is a scalar Mobius transform):
#   c   : shift, H' = H + cI SPD on the sector (c=1.0; verified lambda_min(H')>0).
#   v1  : (e,e,e) Ritz ground vector at PRIMARY resolution N, zero-padded into the
#         AUXILIARY (e,e,e) space at N' >= N (nested cosine modes => exact embedding).
#   A0  = <v1, H' v1>                                     (interval matvec + dot)
#   A2  = <v1, H'^{-1} v1>  via the SPD Goerisch identity, needs NO interval inverse:
#         float-CG solve H'_mid w = v1; r = v1 - H' w (interval residual);
#         A2 = 2<w,v1> - <w,H'w> + <r,H'^{-1}r>,  <r,H'^{-1}r> in [0, ||r||^2 / lam_min],
#         lam_min = c + (rigorous lower bound on mu1(N') = smallest eig of H_{N'}).
#   rho = certified inf(L2) from Stage-A  (rigorous lower bound on lambda_2, and > lambda_1)
#         => rho is a CERTIFIED separator; rho' = rho + c.
#   Mobius (A1=1): A=A0-rho', B=A0-2rho'+rho'^2 A2, nu=A/B,
#         L1_LG = rho' - rho'/(1-nu) - c.   Valid iff inf(B)>0 and sup(nu)<1.
#   inf(L1_LG interval) is the rigorous lower bound on lambda_1.
module LGVerified
using LinearAlgebra, Printf
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
include("assembly_verified.jl")
using .AssemblyVerified
const AV = AssemblyVerified
using .AssemblyVerified.MomentsVerified
using Arpack, IterativeSolvers

const IV = Interval{Float64}

# primary-N eee block linear indices inside the aux-N' eee array (x fastest).
function primary_block_indices(N::Int, Nprime::Int)
    nxs = length(AV.sector_idx(N,:even));  nys = nxs; nzs = nxs
    nxb = length(AV.sector_idx(Nprime,:even)); nyb = nxb; nzb = nxb
    idx = Int[]
    for cz in 1:nzs, cy in 1:nys, cx in 1:nxs
        push!(idx, (cz-1)*nyb*nxb + (cy-1)*nxb + cx)
    end
    return idx, (nxs,nys,nzs), (nxb,nyb,nzb)
end

# interval matvec y = M*x  (M::Matrix{IV}, x::Vector{IV})
function imatvec(M::Matrix{IV}, x::Vector{IV})
    D = length(x); y = Vector{IV}(undef, D)
    @inbounds for i in 1:D
        acc = interval(0.0)
        @simd for j in 1:D
            acc += M[i,j]*x[j]
        end
        y[i] = acc
    end
    return y
end
idot(a::Vector{IV}, b::Vector{IV}) = (s=interval(0.0); @inbounds for i in eachindex(a); s+=a[i]*b[i]; end; s)

function lg_bound(N::Int, Nprime::Int, c::Float64, rho_L2cert::Float64, mu1_lo_Nprime::Float64;
                  nt::Int=48, npanel::Int=96, nodes::Vector{IV}, wts::Vector{IV}, verbose=true)
    @assert Nprime >= N
    T0 = time()
    # --- assemble AUXILIARY eee interval Hamiltonian at N' ---
    ta = @elapsed (P,Kd,Db) = AV.assemble_PK(Nprime,:even,:even,:even; nt=nt, npanel=npanel, nodes=nodes, wts=wts)
    verbose && @printf("  [aux assemble N'=%d] D=%d (%.1fs)\n", Nprime, Db, ta); flush(stdout)
    Hc = AV.form_H(P, Kd; sigma=interval(c))    # H' = H + cI  (interval)
    P = nothing
    Hm = mid.(Hc)                                # float mid for eigensolve + CG

    # --- primary-N Ritz ground vector, extracted from the principal block, zero-padded ---
    pidx, (nxs,nys,nzs), (nxb,nyb,nzb) = primary_block_indices(N, Nprime)
    Hsub = Hm[pidx, pidx]
    es = eigen(Symmetric(Hsub))
    v_small = es.vectors[:,1]
    v = zeros(Float64, Db); v[pidx] .= v_small
    v ./= norm(v)                                # normalized test vector in aux space
    vI = interval.(v)

    # --- A0 = <v, H' v> ---
    Hv = imatvec(Hc, vI)
    A0 = idot(vI, Hv)

    # --- Goerisch: float-CG solve H'_mid w = v ; interval residual & SPD enclosure of A2 ---
    w = zeros(Float64, Db)
    cg!(w, Symmetric(Hm), v; reltol=1e-12, maxiter=8000)
    cg_res = norm(Hm*w .- v)/norm(v)
    wI = interval.(w)
    Hw = imatvec(Hc, wI)
    r = vI .- Hw                                  # interval residual
    wv = idot(wI, vI)
    wHw = idot(wI, Hw)
    r2  = idot(r, r)                              # ||r||^2 enclosure (>=0)
    lam_min_lb = c + mu1_lo_Nprime                # rigorous lower bound on lambda_min(H')
    @assert lam_min_lb > 0 "shift too small: lam_min_lb=$lam_min_lb"
    corr = interval(0.0, sup(r2)/lam_min_lb)      # <r,H'^{-1}r> in [0, ||r||^2/lam_min]
    A2 = interval(2.0)*wv - wHw + corr

    # --- Mobius transform (A1 = 1) ---
    rho_p = interval(rho_L2cert) + interval(c)
    A = A0 - rho_p
    B = A0 - interval(2.0)*rho_p + rho_p*rho_p*A2
    nu = A/B
    lam_hat = rho_p - rho_p/(interval(1.0)-nu)
    L1_LG = lam_hat - interval(c)

    Bpos = inf(B) > 0
    nult1 = sup(nu) < 1
    verbose && @printf("  [LG] A0=[%.12f,%.12f] A2=[%.12f,%.12f]\n", inf(A0),sup(A0),inf(A2),sup(A2))
    verbose && @printf("  [LG] B=[%.4e,%.4e] (pos=%s) nu=[%.5f,%.5f] (<1=%s) cg_res=%.2e\n",
                       inf(B),sup(B),Bpos,inf(nu),sup(nu),nult1,cg_res); flush(stdout)
    verbose && @printf("  [LG] L1_LG=[%.12f,%.12f] width=%.3e\n", inf(L1_LG),sup(L1_LG),diam(L1_LG)); flush(stdout)

    return Dict{String,Any}(
        "N"=>N, "Nprime"=>Nprime, "c"=>c, "D_aux"=>Db,
        "rho_L2cert"=>rho_L2cert, "mu1_lo_Nprime"=>mu1_lo_Nprime, "lam_min_lb"=>lam_min_lb,
        "A0"=>[inf(A0),sup(A0)], "A2"=>[inf(A2),sup(A2)], "A2_corr_hi"=>sup(corr),
        "B"=>[inf(B),sup(B)], "nu"=>[inf(nu),sup(nu)], "B_positive"=>Bpos, "nu_lt_1"=>nult1,
        "cg_residual"=>cg_res, "r2_hi"=>sup(r2),
        "L1_LG"=>[inf(L1_LG),sup(L1_LG)], "L1_LG_lo"=>inf(L1_LG), "L1_LG_width"=>diam(L1_LG),
        "t_assemble"=>ta, "wall_total"=>time()-T0)
end
end # module
