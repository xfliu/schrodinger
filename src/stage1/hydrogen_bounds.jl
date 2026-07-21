#!/usr/bin/env julia
"""
Verified certified lower bounds for the hydrogen-atom ground state — Julia Stage 1.

Operator  H = -Δ - 1/|x|  on the cube  Q = (-6,6)^3  with Neumann BC,
in the (-Δ) convention where the exact ground state is λ₁ = -1/4.

Faithful double-precision port of hydrogen_bounds.py:
  * spectral Galerkin in the L2-normalized Neumann cosine basis
        φ_m(x) = ∏_i c_{m_i}(x_i),  m_i = 0..N,  ν_m = Σ_i (m_i π / 2L)^2
  * Coulomb matrix assembled EXACTLY via the Gaussian (Laplace) representation
        1/r = (2/√π) ∫_0^∞ e^{-t² r²} dt,
    factorized into a Kronecker product of 1D moments
        B_t(n,n') = ∫_{-L}^{L} c_n c_{n'} e^{-t² x²} dx.
    P = -pref Σ_t w_t (B_t ⊗ B_t ⊗ B_t)  → matrix-free matvec.
  * lowest Galerkin eigenvalue μ_{1,N} by preconditioned LOBPCG.
  * certified lower bound L1 from Thm. coulomb-hardy (same formula as Python).
  fixed ρ=6, Z=1, 36/ρ²=1, ε=0.45, σ=4.50 (manuscript).
"""
module HydrogenBounds

using LinearAlgebra
using SpecialFunctions: erfcx
using IterativeSolvers: lobpcg
using Printf
using JSON

const L    = 6.0
const Z    = 1.0
const RHO  = 6.0
const H36  = 36.0 / RHO^2       # = 1.0
const EPS  = 0.45               # manuscript fixed value
const LAM1 = -0.25              # exact hydrogen ground state

# ---- Gauss-Legendre nodes/weights on [-1,1] via Golub-Welsch --------------
function gauss_legendre(n::Int)
    # symmetric tridiagonal Jacobi matrix; β_k = k/√(4k²-1)
    k = 1:(n-1)
    β = k ./ sqrt.(4.0 .* k.^2 .- 1.0)
    J = SymTridiagonal(zeros(n), collect(β))
    vals, vecs = eigen(J)
    w = 2.0 .* (vecs[1, :] .^ 2)
    return vals, w
end

# map GL on [-1,1] to (0,∞) via u=(x+1)/2, t=u/(1-u)
function build_t_grid(nt::Int)
    x, w = gauss_legendre(nt)
    u  = 0.5 .* (x .+ 1); wu = 0.5 .* w
    t  = u ./ (1 .- u)
    wt = wu ./ (1 .- u).^2
    return t, wt
end

# 1D Gaussian moment kernel  g(κ,t) = ∫_{-∞}^{∞}... stable via erfcx (Faddeeva).
# Matches Python g1d: returns array (nκ × nt).
# g = (√π/t)[ e^{-b²} - Re( e^{-a²} e^{-2 i a b} w(i z) ) ],  z=a+ i b, a=tL, b=κ/(2t)
# Using w(iz) = erfcx(z) for real z≥0 -> the whole bracket reduces to a real,
# overflow-free expression:  e^{-b²} - e^{-a²}·Re(e^{-2iab} w(i(a+ib))).
# For real a,b: i z = -b + i a, and w(-b + i a) with Re= -b handled via
# erfcx of complex arg. SpecialFunctions.erfcx accepts Complex.
function g1d(kappa::AbstractVector, t::AbstractVector, L::Float64)
    nκ = length(kappa); nt = length(t)
    G = Array{Float64}(undef, nκ, nt)
    @inbounds for j in 1:nt
        tj = t[j]
        for i in 1:nκ
            a = tj * L
            b = kappa[i] / (2tj)
            z = a + im*b                 # matches Python zz = a + i b
            # Python: (√π/t)*( e^{-b²} - Re( e^{-a²} e^{-2 i a b} wofz(i z) ) )
            # wofz(iz) = erfcx(z)  (since w(iy)=erfcx(y); here arg is i z)
            wv = erfcx(z)                # = wofz(i z)
            term = exp(-a^2) * exp(-2im*a*b) * wv
            G[i, j] = (sqrt(pi)/tj) * (exp(-b^2) - real(term))
        end
    end
    return G
end

# 1D moment tensor B[t,n,m], n,m=0..N  (returns Array nt × (N+1) × (N+1))
function Bt_tensor(N::Int, L::Float64, t::AbstractVector)
    fmax = 2N
    freqs = collect(0:fmax) .* (pi/(2L))
    GC = g1d(freqs, t, L)                       # (fmax+1) × nt
    coef = cos.(collect(0:fmax) .* (pi/2))      # = 0 for odd, ±1 for even
    Nn = [ (n == 0) ? 1/sqrt(2L) : 1/sqrt(L) for n in 0:N ]
    nt = length(t)
    B = zeros(Float64, nt, N+1, N+1)
    @inbounds for n in 0:N, m in 0:N
        fm = abs(n-m); fp = n+m
        for k in 1:nt
            B[k, n+1, m+1] = Nn[n+1]*Nn[m+1]*0.5*(coef[fm+1]*GC[fm+1,k] + coef[fp+1]*GC[fp+1,k])
        end
    end
    return B
end

# matrix-free block matmat: X is D×kk (D=(N+1)^3). Mirrors Python einsum path.
# Y = -pref Σ_t w_t (B_t ⊗ B_t ⊗ B_t) X  +  Kdiag .* X
function make_block_matmat(N::Int, L::Float64, t::AbstractVector, wt::AbstractVector; sigma=0.0)
    B = Bt_tensor(N, L, t)                       # nt × n1 × n1
    n1 = N+1; nt = length(t)
    pref = 2/sqrt(pi)
    Wt = wt .* (-pref)
    nu = (pi/(2L))^2
    m = collect(0:N)
    Kdiag = zeros(Float64, n1^3)
    idx = 1
    @inbounds for a in 0:N, b in 0:N, c in 0:N
        Kdiag[idx] = nu*(a^2+b^2+c^2) + sigma
        idx += 1
    end
    # per-t 1D operators as matrices for fast mode-multiplication
    Bmats = [ Matrix{Float64}(B[k, :, :]) for k in 1:nt ]  # each n1×n1, symmetric

    # apply (Bk ⊗ Bk ⊗ Bk) to a single vector reshaped (n1,n1,n1)
    # via three mode multiplications. Index order matches Python C=reshape(kk,n1,n1,n1)
    # with axes (i,j,l): mode-1 over i, mode-2 over j, mode-3 over l.
    function matmat(X::AbstractMatrix)
        D, kk = size(X)
        R = zeros(Float64, D, kk)
        tmp = Array{Float64}(undef, n1, n1, n1)
        for col in 1:kk
            C = reshape(view(X, :, col), n1, n1, n1)
            acc = zeros(Float64, n1, n1, n1)
            for k in 1:nt
                Bk = Bmats[k]
                # mode-1: contract Bk over first index
                Y1 = reshape(Bk * reshape(C, n1, n1*n1), n1, n1, n1)
                # mode-2: over second index
                Y2 = permutedims(Y1, (2,1,3))
                Y2 = reshape(Bk * reshape(Y2, n1, n1*n1), n1, n1, n1)
                Y2 = permutedims(Y2, (2,1,3))
                # mode-3: over third index
                Y3 = permutedims(Y2, (3,1,2))
                Y3 = reshape(Bk * reshape(Y3, n1, n1*n1), n1, n1, n1)
                Y3 = permutedims(Y3, (2,3,1))
                acc .+= Wt[k] .* Y3
            end
            R[:, col] = reshape(acc, D) .+ Kdiag .* view(X, :, col)
        end
        return R
    end
    matvec(x::AbstractVector) = vec(matmat(reshape(x, :, 1)))
    return matmat, matvec, Kdiag
end

# A LinearMap-like struct implementing mul! for IterativeSolvers.lobpcg
struct BlockOp
    matmat::Function
    D::Int
end
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

# Diagonal preconditioner P⁻¹ ≈ 1/(Kdiag+5)
struct DiagPrec
    d::Vector{Float64}
end
function mul!(Y::AbstractVecOrMat, P::DiagPrec, X::AbstractVecOrMat)
    Y .= P.d .* X
    return Y
end
Base.size(P::DiagPrec) = (length(P.d), length(P.d))
Base.size(P::DiagPrec, i::Int) = length(P.d)
Base.eltype(::DiagPrec) = Float64
LinearAlgebra.issymmetric(::DiagPrec) = true
ldiv!(Y, P::DiagPrec, X) = (Y .= P.d .* X; Y)

function galerkin_mu(N::Int; nt::Int=48, k::Int=2, tol::Float64=1e-8, maxiter::Int=400)
    t, wt = build_t_grid(nt)
    matmat, matvec, Kd = make_block_matmat(N, L, t, wt)
    D = (N+1)^3
    A = BlockOp(matmat, D)
    Pinv = DiagPrec(1.0 ./ (Kd .+ 5.0))
    r = lobpcg(A, false, k; P=Pinv, tol=tol, maxiter=maxiter)
    vals = sort(r.λ)
    return vals[1:k]
end

function certified_L1(N::Int, mu1N::Float64; eps::Float64=EPS)
    sigma = (eps/8)*H36 + 2*Z^2/eps
    ms = mu1N + sigma
    etaV = Z*sqrt(8/(1-eps) + H36/ms)
    g = ms^(-0.5)*etaV
    nu_star = (N+1)^2*pi^2/(4*L^2)
    Chat2 = (1+g^2)/((1-eps)*nu_star)
    L1 = ms/(1+Chat2*ms) - sigma
    return Dict("N"=>N, "mu1N"=>mu1N, "sigma"=>sigma, "g"=>g,
                "nu_star"=>nu_star, "Chat2"=>Chat2,
                "L1"=>L1, "gap"=>LAM1-L1)
end

function main(Ns::Vector{Int})
    nt = 48
    results = Dict[]
    for N in Ns
        t0 = time()
        ev = galerkin_mu(N; nt=nt)
        mu1 = ev[1]; mu2 = ev[2]
        rec = certified_L1(N, mu1)
        rec["mu2N"] = mu2; rec["dofs"] = (N+1)^3; rec["seconds"] = time()-t0
        push!(results, rec)
        @printf("N=%4d D=%9d mu1=%.8f sigma=%.4f nu*=%.2f L1=%.6f gap=%.6f (%.1fs)\n",
                N, (N+1)^3, mu1, rec["sigma"], rec["nu_star"], rec["L1"], rec["gap"], rec["seconds"])
        flush(stdout)
    end
    for i in 2:length(results)
        g0 = results[i-1]["gap"]; g1 = results[i]["gap"]
        N0 = results[i-1]["N"]; N1 = results[i]["N"]
        if g0 > 0 && g1 > 0
            results[i]["EOC"] = log(g0/g1)/log(N1/N0)
        end
    end
    open("hydrogen_bounds_julia.json", "w") do f
        JSON.print(f, results, 2)
    end
    println("WROTE hydrogen_bounds_julia.json")
    return results
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Ns = length(ARGS) > 0 ? parse.(Int, ARGS) : [16,24,32,48,64,96,128]
    HydrogenBounds.main(Ns)
end
