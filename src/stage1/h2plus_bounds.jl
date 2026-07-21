#!/usr/bin/env julia
"""
H2+ molecular ion — two-centre Coulomb, Julia Stage 1 (double precision).

Operator  H = -Δ - 1/|x-a1| - 1/|x-a2|,  a_{1,2}=(∓2,0,0),  on the enlarged
anisotropic box  Q' = [-10,10]×[-8,8]^2  (Neumann), -Δ convention.
References  λ1 = -0.551317,  λ2 = -0.333768 ;  U1 ≈ -0.545.

Spectral Galerkin in the L2-normalized Neumann cosine basis, per-axis:
    φ_m(x)=∏_i c_{m_i}(x_i;L_i),  L=(10,8,8),  ν_m=Σ_i (m_i π/2L_i)^2.
Coulomb assembled EXACTLY via Gaussian (Laplace) representation
    1/|x-a| = (2/√π)∫_0^∞ e^{-t²|x-a|²} dt,  e^{-t²|x-a|²}=∏_i e^{-t²(x_i-a_i)²},
factorized into a Kronecker product of SHIFTED 1D moments
    B_t^{(s)}(n,m)=∫_{-L}^{L} c_n c_m e^{-t²(x-s)²} dx.
Two nuclei share the y,z factors (a=(∓2,0,0), unshifted in y,z); only the
x-factor is shifted, by s=∓2. So  Bx_t = Bx_t^{(-2)} + Bx_t^{(+2)}.
Matrix-free matvec via three mode multiplications (distinct op per axis).
Certified bounds L1,L2 from thm:coulomb-hardy, Z_tot=2, ρ=8, ε optimised per N.
"""
module H2plusBounds

using LinearAlgebra
using SpecialFunctions: erfcx
using IterativeSolvers: lobpcg
using Printf
using JSON

const LX, LY, LZ = 10.0, 8.0, 8.0
const AX   = 2.0                    # nuclei at (∓AX,0,0)
const ZTOT = 2.0
const RHO  = 8.0                    # min dist(nucleus, ∂Q'): (∓2,0,0)->y,z faces at ±8
const H36  = 36.0 / RHO^2           # = 0.5625
const LAM1 = -0.551317
const LAM2 = -0.333768

# ---- Gauss-Legendre via Golub-Welsch --------------------------------------
function gauss_legendre(n::Int)
    k = 1:(n-1)
    β = k ./ sqrt.(4.0 .* k.^2 .- 1.0)
    J = SymTridiagonal(zeros(n), collect(β))
    vals, vecs = eigen(J)
    return vals, 2.0 .* (vecs[1, :] .^ 2)
end
function build_t_grid(nt::Int)
    x, w = gauss_legendre(nt)
    u  = 0.5 .* (x .+ 1); wu = 0.5 .* w
    return u ./ (1 .- u), wu ./ (1 .- u).^2
end

# ---- SHIFTED 1D cosine moment, Faddeeva-stable (reflection-bounded) --------
# B[k,n,m] = ∫_{-L}^{L} c_n(x;L) c_m(x;L) e^{-t²(x-s)²} dx,  shift s.
# c_n = Nn cos(κ_n (x+L)),  κ_n=nπ/2L,  Nn=1/√(2L) (n=0) else 1/√L.
# Per frequency κ:  G_κ = Re[ e^{iκ(L+s)} (√π/2t)(Bnd(α)-Bnd(β)) ],  α=-L-s, β=L-s,
# Bnd(y0)= (y0≥0) ? erfcx(z)·base  :  2e^{-κ²/4t²} - erfcx(-z)·base,
#   z = t·y0 - iκ/(2t),  base = e^{-t²y0² + iκ y0}.   (wofz(iz)=erfcx(z))
function shifted_moment(N::Int, L::Float64, s::Float64, t::AbstractVector)
    n1 = N+1; nt = length(t); fmax = 2N
    κ  = collect(0:fmax) .* (pi/(2L))
    Nn = [ (n==0) ? 1/sqrt(2L) : 1/sqrt(L) for n in 0:N ]
    α  = -L - s; β = L - s
    B  = zeros(Float64, nt, n1, n1)
    G  = Vector{Float64}(undef, fmax+1)
    @inbounds for k in 1:nt
        tt = t[k]
        for (fi, kap) in enumerate(κ)
            # Bnd(α)
            za = tt*α - im*kap/(2tt); basea = exp(-tt^2*α^2 + im*kap*α)
            Ba = α >= 0 ? erfcx(za)*basea : 2*exp(-kap^2/(4tt^2)) - erfcx(-za)*basea
            # Bnd(β)
            zb = tt*β - im*kap/(2tt); baseb = exp(-tt^2*β^2 + im*kap*β)
            Bb = β >= 0 ? erfcx(zb)*baseb : 2*exp(-kap^2/(4tt^2)) - erfcx(-zb)*baseb
            integ = (sqrt(pi)/(2tt))*(Ba - Bb)
            G[fi] = real(exp(im*kap*(L+s))*integ)
        end
        for n in 0:N, m in 0:N
            fm = abs(n-m); fp = n+m
            B[k, n+1, m+1] = 0.5*Nn[n+1]*Nn[m+1]*(G[fm+1] + G[fp+1])
        end
    end
    return B
end

# ---- matrix-free matmat with DISTINCT operator per axis --------------------
# Ax = -pref Σ_t w_t (Bx_t ⊗ By_t ⊗ Bz_t) x  +  Kdiag .* x
function make_matmat(N::Int, t::AbstractVector, wt::AbstractVector)
    n1 = N+1; nt = length(t); pref = 2/sqrt(pi); Wt = wt .* (-pref)
    Bxm = shifted_moment(N, LX, -AX, t); Bxp = shifted_moment(N, LX, +AX, t)
    Byt = shifted_moment(N, LY, 0.0, t); Bzt = shifted_moment(N, LZ, 0.0, t)
    Bx = [ Matrix{Float64}(Bxm[k,:,:] .+ Bxp[k,:,:]) for k in 1:nt ]
    By = [ Matrix{Float64}(Byt[k,:,:]) for k in 1:nt ]
    Bz = [ Matrix{Float64}(Bzt[k,:,:]) for k in 1:nt ]
    νx = (pi/(2LX))^2; νy = (pi/(2LY))^2; νz = (pi/(2LZ))^2
    Kdiag = zeros(Float64, n1^3); idx = 1
    # Julia reshape is COLUMN-MAJOR: axis-1 (i, the Bx axis) varies FASTEST.
    # So the innermost loop variable must carry νx to stay consistent with matvec.
    @inbounds for c in 0:N, b in 0:N, a in 0:N
        Kdiag[idx] = νx*a^2 + νy*b^2 + νz*c^2; idx += 1
    end
    function matmat(X::AbstractMatrix)
        D, kk = size(X); R = zeros(Float64, D, kk)
        for col in 1:kk
            C = reshape(view(X, :, col), n1, n1, n1)
            acc = zeros(Float64, n1, n1, n1)
            for k in 1:nt
                # mode-1 (Bx over axis i)
                Y1 = reshape(Bx[k] * reshape(C, n1, n1*n1), n1, n1, n1)
                # mode-2 (By over axis j)
                Y2 = permutedims(Y1, (2,1,3))
                Y2 = reshape(By[k] * reshape(Y2, n1, n1*n1), n1, n1, n1)
                Y2 = permutedims(Y2, (2,1,3))
                # mode-3 (Bz over axis l)
                Y3 = permutedims(Y2, (3,1,2))
                Y3 = reshape(Bz[k] * reshape(Y3, n1, n1*n1), n1, n1, n1)
                Y3 = permutedims(Y3, (2,3,1))
                acc .+= Wt[k] .* Y3
            end
            R[:, col] = reshape(acc, D) .+ Kdiag .* view(X, :, col)
        end
        return R
    end
    return matmat, Kdiag
end

struct BlockOp; matmat::Function; D::Int; end
Base.size(A::BlockOp) = (A.D, A.D)
Base.size(A::BlockOp, i::Int) = A.D
Base.eltype(::BlockOp) = Float64
LinearAlgebra.issymmetric(::BlockOp) = true
import LinearAlgebra: mul!, ldiv!
function mul!(Y::AbstractVecOrMat, A::BlockOp, X::AbstractVecOrMat)
    if X isa AbstractVector
        Y .= vec(A.matmat(reshape(X, :, 1)))
    else
        Y .= A.matmat(X)
    end
    return Y
end
Base.:*(A::BlockOp, X::AbstractVecOrMat) = A.matmat(X isa AbstractVector ? reshape(X,:,1) : X)

struct DiagPrec; d::Vector{Float64}; end
function mul!(Y::AbstractVecOrMat, P::DiagPrec, X::AbstractVecOrMat); Y .= P.d .* X; Y; end
Base.size(P::DiagPrec) = (length(P.d), length(P.d))
Base.size(P::DiagPrec, i::Int) = length(P.d)
Base.eltype(::DiagPrec) = Float64
LinearAlgebra.issymmetric(::DiagPrec) = true
ldiv!(Y, P::DiagPrec, X) = (Y .= P.d .* X; Y)

function galerkin_mu(N::Int; nt::Int=48, k::Int=3, tol::Float64=1e-8, maxiter::Int=400)
    t, wt = build_t_grid(nt)
    matmat, Kd = make_matmat(N, t, wt)
    D = (N+1)^3
    A = BlockOp(matmat, D)
    Pinv = DiagPrec(1.0 ./ (Kd .+ 5.0))
    r = lobpcg(A, false, k; P=Pinv, tol=tol, maxiter=maxiter)
    return sort(r.λ)[1:k]
end

# ---- certified Hardy bounds for k=1,2, ε optimised to maximise L2 ----------
function Cconst(delta::Float64, r1::Float64)
    # Localized 3D Hardy: chi=1 on B(a,r1), supp in B(a,RHO), |grad chi|<=1/(RHO-r1)
    #   int |v|^2/r^2 <= Cg ||grad v||^2 + CL ||v||^2 (per centre)
    # Young split (1+delta): Cg = 4(1+delta);  CL = 4(1+1/delta)/(RHO-r1)^2 + 1/r1^2
    # Manuscript special case delta=1, r1=RHO/2 -> Cg=8, CL=36/RHO^2.
    Cg = 4.0*(1.0+delta)
    CL = 4.0*(1.0+1.0/delta)/(RHO-r1)^2 + 1.0/r1^2
    return Cg, CL
end

function hardy_L(mu1N::Float64, mukN::Float64, N::Int, eps::Float64,
                 delta::Float64=1.0, r1::Float64=RHO/2)
    Cg, CL = Cconst(delta, r1)
    sigma = eps*CL/Cg + ZTOT^2*Cg/(4*eps)
    ms1 = mu1N + sigma
    etaV = ZTOT*sqrt(Cg/(1-eps) + CL/ms1)
    g = ms1^(-0.5)*etaV
    nu = ((N+1)*pi/(2LX))^2
    Chat2 = (1+g^2)/((1-eps)*nu)
    msk = mukN + sigma
    return msk/(1+Chat2*msk) - sigma, sigma, g, nu, Chat2
end
function optimize_eps(mu1N, mukN, N; ne=400, nd=120, nr=120)
    # Joint optimization over (delta, r1, eps) maximizing the certified L.
    # Coarse grid then local golden refine on eps at the best (delta,r1).
    best=-Inf; be=0.45; bd=1.0; br=RHO/2
    for delta in exp.(range(log(5e-3), log(3.0); length=nd))
        for r1 in range(0.05, RHO-0.05; length=nr)
            Cg, CL = Cconst(delta, r1)
            for e in range(1e-3, 0.999; length=ne)
                sigma = e*CL/Cg + ZTOT^2*Cg/(4*e)
                if mu1N + sigma <= 0; continue; end
                L,_ = hardy_L(mu1N, mukN, N, e, delta, r1)
                if L>best; best=L; be=e; bd=delta; br=r1; end
            end
        end
    end
    # refine eps on a fine local grid around be
    lo=max(1e-3, be-0.02); hi=min(0.999, be+0.02)
    for e in range(lo, hi; length=4000)
        L,_ = hardy_L(mu1N, mukN, N, e, bd, br)
        if L>best; best=L; be=e; end
    end
    return be, best, bd, br
end

function main(Ns::Vector{Int}; nt::Int=48)
    results = Dict[]
    for N in Ns
        t0 = time()
        ev = galerkin_mu(N; nt=nt, k=3)
        mu1, mu2 = ev[1], ev[2]
        e1, L1, d1, r1a = optimize_eps(mu1, mu1, N)
        e2, L2, d2, r2a = optimize_eps(mu1, mu2, N)
        nu = ((N+1)*pi/(2LX))^2
        rec = Dict("N"=>N, "dofs"=>(N+1)^3, "mu1N"=>mu1, "mu2N"=>mu2,
                   "mu3N"=>ev[3], "eps1"=>e1, "eps2"=>e2, "L1"=>L1, "L2"=>L2,
                   "nu_star"=>nu, "gap2"=>LAM2-L2, "sep_L2_gt_U1"=>(L2 > -0.545),
                   "delta2"=>d2, "r1_2"=>r2a, "delta1"=>d1, "r1_1"=>r1a,
                   "seconds"=>time()-t0)
        push!(results, rec)
        @printf("N=%4d D=%10d mu1=%.8f mu2=%.8f nu*=%.1f L1=%.4f L2=%.4f gap2=%.4f sep=%s (%.1fs)\n",
                N, (N+1)^3, mu1, mu2, nu, L1, L2, LAM2-L2, (L2 > -0.545), rec["seconds"])
        flush(stdout)
    end
    for i in 2:length(results)
        g0=results[i-1]["gap2"]; g1=results[i]["gap2"]
        N0=results[i-1]["N"]; N1=results[i]["N"]
        if g0>0 && g1>0; results[i]["EOC"]=log(g0/g1)/log(N1/N0); end
    end
    open("h2plus_bounds_julia.json","w") do f; JSON.print(f, results, 2); end
    println("WROTE h2plus_bounds_julia.json")
    return results
end


# ---- parity-sector matvec: restrict each axis to even (0-based even idx) or odd
function make_matmat_sector(N::Int, t::AbstractVector, wt::AbstractVector,
                            px::Symbol, py::Symbol, pz::Symbol)
    n1 = N+1; nt = length(t); pref = 2/sqrt(pi); Wt = wt .* (-pref)
    Bxm = shifted_moment(N, LX, -AX, t); Bxp = shifted_moment(N, LX, +AX, t)
    Byt = shifted_moment(N, LY, 0.0, t); Bzt = shifted_moment(N, LZ, 0.0, t)
    sel(p) = p === :even ? collect(1:2:n1) : collect(2:2:n1)   # 1-based indices
    Ix, Iy, Iz = sel(px), sel(py), sel(pz)
    Bx = [ Matrix{Float64}((Bxm[k,:,:] .+ Bxp[k,:,:])[Ix, Ix]) for k in 1:nt ]
    By = [ Matrix{Float64}(Byt[k,:,:][Iy, Iy]) for k in 1:nt ]
    Bz = [ Matrix{Float64}(Bzt[k,:,:][Iz, Iz]) for k in 1:nt ]
    νx = (pi/(2LX))^2; νy = (pi/(2LY))^2; νz = (pi/(2LZ))^2
    nx, ny, nz = length(Ix), length(Iy), length(Iz)
    mx = Ix .- 1; my = Iy .- 1; mz = Iz .- 1     # 0-based mode numbers
    Kdiag = zeros(Float64, nx*ny*nz); idx = 1
    @inbounds for c in 1:nz, b in 1:ny, a in 1:nx
        Kdiag[idx] = νx*mx[a]^2 + νy*my[b]^2 + νz*mz[c]^2; idx += 1
    end
    function matmat(X::AbstractMatrix)
        D, kk = size(X); R = zeros(Float64, D, kk)
        for col in 1:kk
            C = reshape(view(X, :, col), nx, ny, nz)
            acc = zeros(Float64, nx, ny, nz)
            for k in 1:nt
                Y1 = reshape(Bx[k] * reshape(C, nx, ny*nz), nx, ny, nz)
                Y2 = permutedims(Y1, (2,1,3))
                Y2 = reshape(By[k] * reshape(Y2, ny, nx*nz), ny, nx, nz)
                Y2 = permutedims(Y2, (2,1,3))
                Y3 = permutedims(Y2, (3,1,2))
                Y3 = reshape(Bz[k] * reshape(Y3, nz, nx*ny), nz, nx, ny)
                Y3 = permutedims(Y3, (2,3,1))
                acc .+= Wt[k] .* Y3
            end
            R[:, col] = reshape(acc, D) .+ Kdiag .* view(X, :, col)
        end
        return R
    end
    return matmat, Kdiag, nx*ny*nz
end

function galerkin_sector(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                         nt::Int=48, k::Int=1, tol::Float64=1e-8, maxiter::Int=400)
    t, wt = build_t_grid(nt)
    matmat, Kd, D = make_matmat_sector(N, t, wt, px, py, pz)
    A = BlockOp(matmat, D)
    Pinv = DiagPrec(1.0 ./ (Kd .+ 5.0))
    r = lobpcg(A, false, k; P=Pinv, tol=tol, maxiter=maxiter)
    return sort(r.λ)[1:k], D
end

# μ1 = ground of (even,even,even);  μ2 = ground of (odd,even,even)
function main_sector(Ns::Vector{Int}; nt::Int=48, tol::Float64=1e-8, maxiter::Int=800)
    results = Dict[]
    for N in Ns
        t0 = time()
        ev1, D1 = galerkin_sector(N, :even,:even,:even; nt=nt, k=1, tol=tol, maxiter=maxiter)
        ev2, D2 = galerkin_sector(N, :odd, :even,:even; nt=nt, k=1, tol=tol, maxiter=maxiter)
        mu1, mu2 = ev1[1], ev2[1]
        e1, L1, d1, r1a = optimize_eps(mu1, mu1, N)
        e2, L2, d2, r2a = optimize_eps(mu1, mu2, N)
        nu = ((N+1)*pi/(2LX))^2
        rec = Dict("N"=>N, "dofs_full"=>(N+1)^3, "dofs_eee"=>D1, "dofs_oee"=>D2,
                   "mu1N"=>mu1, "mu2N"=>mu2, "eps1"=>e1, "eps2"=>e2,
                   "L1"=>L1, "L2"=>L2, "nu_star"=>nu, "gap2"=>LAM2-L2,
                   "sep_L2_gt_U1"=>(L2 > -0.545),
                   "delta2"=>d2, "r1_2"=>r2a, "delta1"=>d1, "r1_1"=>r1a,
                   "seconds"=>time()-t0)
        push!(results, rec)
        @printf("N=%4d Deee=%9d Doee=%9d mu1=%.8f mu2=%.8f nu*=%.1f L1=%.4f L2=%.4f gap2=%.4f sep=%s (%.1fs)\n",
                N, D1, D2, mu1, mu2, nu, L1, L2, LAM2-L2, (L2 > -0.545), rec["seconds"])
        flush(stdout)
    end
    for i in 2:length(results)
        g0=results[i-1]["gap2"]; g1=results[i]["gap2"]; N0=results[i-1]["N"]; N1=results[i]["N"]
        if g0>0 && g1>0; results[i]["EOC"]=log(g0/g1)/log(N1/N0); end
    end
    open("h2plus_bounds_sector.json","w") do f; JSON.print(f, results, 2); end
    println("WROTE h2plus_bounds_sector.json")
    return results
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    args = ARGS
    if length(args) > 0 && args[1] == "sector"
        Ns = length(args) > 1 ? parse.(Int, args[2:end]) : [24,32,48,64]
        H2plusBounds.main_sector(Ns)
    else
        Ns = length(args) > 0 ? parse.(Int, args) : [12,16,24,32,48,64,96,128]
        H2plusBounds.main(Ns)
    end
end
