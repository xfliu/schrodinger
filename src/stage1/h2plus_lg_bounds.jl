#!/usr/bin/env julia
"""
H2+ molecular ion — single-test-vector Lehmann–Goerisch (LG) lower bound for
lambda_1, built on top of the Neumann cosine sector-Galerkin machinery in
h2plus_bounds.jl.

Setting (unchanged from h2plus_bounds.jl):
    H = -Δ - 1/|x-a1| - 1/|x-a2|,  a_{1,2}=(∓2,0,0),  on Q'=[-10,10]×[-8,8]^2
    (Neumann cosine basis, per-axis modes 0..N).  Ground state lambda_1 lives
    in the (even,even,even) parity sector; the first excited state lambda_2
    lives in the (odd,even,even) sector (h2plus_bounds.jl notation).

Method (single test vector -> LG reduces to a linear solve, no saddle-point
mixed-FEM system needed):

  1. SHIFT.  Choose c > -mu1(N) so that H' = H + c is SPD on the sector space
     at every resolution used below (verified numerically).

  2. TEST VECTOR.  v1 = normalized (even,even,even) Ritz ground-state
     eigenvector at PRIMARY resolution N (from lobpcg).

  3. AUXILIARY SPACE.  Embed v1, by zero-padding in mode-number, into the
     richer (even,even,even) sector space at AUXILIARY resolution N' >= N.
     Because the cosine-Galerkin matrix elements B_t^{(s)}(n,m) do not depend
     on the truncation level (only on the mode pair (n,m) itself), H_N is
     EXACTLY the leading principal sub-block of H_{N'}; hence
         A0 := <v1, (H_{N'}+c) v1>_{N'} = mu1(N) + c        (no extra call needed)
     This is the discrete analogue of assumption A3 in implementation_LG_method.html
     applied with T = identity embedding (b_G = the H'-inner product itself);
     for a single test vector there is no gradient-recovery / mixed-FEM system —
     A4 is simply "solve H'_{N'} w1 = v1 in the auxiliary space."

  4. GOERISCH STEP (A4).  Solve  (H_{N'}+c) w1 = v1_pad  via CG (SPD after
     shift; matrix-free matvec).  A2 := <v1_pad, w1>.

  5. RHO.  rho = mu2(N) (the (odd,even,even) Ritz value, a numerically stable
     but NOT rigorously certified upper bound for lambda_2 — flagged
     explicitly in the memo).  rho' = rho + c.

  6. MOBIUS TRANSFORM (n=1, closed form).  With A1=1:
         A  = A0 - rho'
         B  = A0 - 2 rho' + rho'^2 A2
         nu = A/B
         lambda_hat_1_LG = rho' - rho'/(1-nu)
         lambda_1_LG     = lambda_hat_1_LG - c
     Requires B>0 (checked) and nu<1 (checked) for the bound to be valid.

Because rho uses mu2(N) rather than a certified lower bound for lambda_2,
the resulting L1_LG is a NUMERICALLY VALIDATED bound (checked against the
oracle choice rho=lam2 in a synthetic test, see h2plus_lg_memo.md), not a
fully certified one in the strict sense of implementation_LG_method.html.
"""
module H2plusLG

using LinearAlgebra
using SpecialFunctions: erfcx
using IterativeSolvers: lobpcg, cg, cg!
using Printf
using JSON

const LX, LY, LZ = 10.0, 8.0, 8.0
const AX   = 2.0
const ZTOT = 2.0

# ---- Gauss-Legendre / shifted cosine moments (identical to h2plus_bounds.jl) ----
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
            za = tt*α - im*kap/(2tt); basea = exp(-tt^2*α^2 + im*kap*α)
            Ba = α >= 0 ? erfcx(za)*basea : 2*exp(-kap^2/(4tt^2)) - erfcx(-za)*basea
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

# ---- parity-sector matrix-free operator (identical structure to h2plus_bounds.jl) ----
function make_matmat_sector(N::Int, t::AbstractVector, wt::AbstractVector,
                            px::Symbol, py::Symbol, pz::Symbol)
    n1 = N+1; nt = length(t); pref = 2/sqrt(pi); Wt = wt .* (-pref)
    Bxm = shifted_moment(N, LX, -AX, t); Bxp = shifted_moment(N, LX, +AX, t)
    Byt = shifted_moment(N, LY, 0.0, t); Bzt = shifted_moment(N, LZ, 0.0, t)
    sel(p) = p === :even ? collect(1:2:n1) : collect(2:2:n1)
    Ix, Iy, Iz = sel(px), sel(py), sel(pz)
    Bx = [ Matrix{Float64}((Bxm[k,:,:] .+ Bxp[k,:,:])[Ix, Ix]) for k in 1:nt ]
    By = [ Matrix{Float64}(Byt[k,:,:][Iy, Iy]) for k in 1:nt ]
    Bz = [ Matrix{Float64}(Bzt[k,:,:][Iz, Iz]) for k in 1:nt ]
    νx = (pi/(2LX))^2; νy = (pi/(2LY))^2; νz = (pi/(2LZ))^2
    nx, ny, nz = length(Ix), length(Iy), length(Iz)
    mx = Ix .- 1; my = Iy .- 1; mz = Iz .- 1
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
    return matmat, Kdiag, nx, ny, nz
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

# ---- shifted operator  H' = H + c*I,  matrix-free, SPD for c large enough ----
struct ShiftedOp; base::BlockOp; c::Float64; end
Base.size(A::ShiftedOp) = size(A.base)
Base.size(A::ShiftedOp, i::Int) = size(A.base, i)
Base.eltype(::ShiftedOp) = Float64
LinearAlgebra.issymmetric(::ShiftedOp) = true
function mul!(Y::AbstractVecOrMat, A::ShiftedOp, X::AbstractVecOrMat)
    mul!(Y, A.base, X)
    Y .+= A.c .* X
    return Y
end
Base.:*(A::ShiftedOp, x::AbstractVector) = A.base*x .+ A.c.*x

struct DiagPrec; d::Vector{Float64}; end
function mul!(Y::AbstractVecOrMat, P::DiagPrec, X::AbstractVecOrMat); Y .= P.d .* X; Y; end
Base.size(P::DiagPrec) = (length(P.d), length(P.d))
Base.size(P::DiagPrec, i::Int) = length(P.d)
Base.eltype(::DiagPrec) = Float64
LinearAlgebra.issymmetric(::DiagPrec) = true
ldiv!(Y, P::DiagPrec, X) = (Y .= P.d .* X; Y)

# ---- ground-state EIGENVECTOR (not just eigenvalue) at resolution N, given sector ----
function ground_state_sector(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                              nt::Int=48, tol::Float64=1e-10, maxiter::Int=800)
    t, wt = build_t_grid(nt)
    matmat, Kd, nx, ny, nz = make_matmat_sector(N, t, wt, px, py, pz)
    D = nx*ny*nz
    A = BlockOp(matmat, D)
    Pinv = DiagPrec(1.0 ./ (Kd .+ 5.0))
    r = lobpcg(A, false, 1; P=Pinv, tol=tol, maxiter=maxiter)
    lam = r.λ[1]
    v = r.X[:, 1]
    v ./= norm(v)
    return lam, v, nx, ny, nz, Kd
end

# ---- zero-pad a (nx,ny,nz) sector coefficient array into a larger (nx',ny',nz') one ----
# Valid because mode lists {0,2,4,...} (even) or {1,3,5,...} (odd) are NESTED with
# increasing resolution: the first nx (resp ny,nz) entries of the big-resolution
# mode list are IDENTICAL to the full small-resolution mode list.
function pad_sector_vector(v_small::Vector{Float64}, nx_s::Int, ny_s::Int, nz_s::Int,
                            nx_b::Int, ny_b::Int, nz_b::Int)
    @assert nx_b >= nx_s && ny_b >= ny_s && nz_b >= nz_s
    C_s = reshape(v_small, nx_s, ny_s, nz_s)
    C_b = zeros(Float64, nx_b, ny_b, nz_b)
    C_b[1:nx_s, 1:ny_s, 1:nz_s] .= C_s
    return reshape(C_b, nx_b*ny_b*nz_b)
end

# ---- single-test-vector LG lower bound for lambda_1 ----
# N       : primary resolution (defines v1 and mu1)
# Nprime  : auxiliary resolution (>= N) used to resolve the Goerisch solve
# c       : shift (must exceed -mu1(N); checked)
# rho_source: :ritz uses mu2(N) (odd,even,even Ritz value); :oracle allows passing rho directly
function lg_bound(N::Int, Nprime::Int, c::Float64; nt::Int=48,
                   cg_tol::Float64=1e-11, cg_maxiter::Int=5000,
                   rho_override::Union{Nothing,Float64}=nothing)
    @assert Nprime >= N
    t0 = time()
    # 1. primary ground state (even,even,even) at resolution N
    mu1, v1, nx_s, ny_s, nz_s, _ = ground_state_sector(N, :even, :even, :even; nt=nt)
    # 2. mu2 from (odd,even,even) at the SAME resolution N (Ritz upper bound for lambda_2)
    mu2, _, _, _, _, _ = ground_state_sector(N, :odd, :even, :even; nt=nt)
    rho = rho_override === nothing ? mu2 : rho_override

    @assert c > -mu1 "shift c=$c not large enough: need c > -mu1(N)=$(-mu1)"

    # 3. auxiliary space operator at resolution N'
    t, wt = build_t_grid(nt)
    matmat_b, Kd_b, nx_b, ny_b, nz_b = make_matmat_sector(Nprime, t, wt, :even, :even, :even)
    Db = nx_b*ny_b*nz_b
    Hbase = BlockOp(matmat_b, Db)
    Hc = ShiftedOp(Hbase, c)

    v1_pad = pad_sector_vector(v1, nx_s, ny_s, nz_s, nx_b, ny_b, nz_b)

    # A0 = <v1_pad, Hc v1_pad> ; should equal mu1+c exactly (principal-submatrix identity) —
    # compute both ways as a consistency check.
    A0_direct = mu1 + c
    Hv1 = Hc * v1_pad
    A0_check = dot(v1_pad, Hv1)

    # 4. Goerisch solve:  Hc * w1 = v1_pad   (SPD system, CG with diagonal preconditioner)
    Pinv = DiagPrec(1.0 ./ (Kd_b .+ c))
    w1 = zeros(Float64, Db)
    hist = cg!(w1, Hc, v1_pad; Pl=Pinv, reltol=cg_tol, maxiter=cg_maxiter, log=true)
    w1_sol, ch = hist
    resid = norm(Hc*w1_sol .- v1_pad) / norm(v1_pad)
    A2 = dot(v1_pad, w1_sol)

    # 5/6. Mobius transform (A1 = 1 since v1 normalized)
    A1 = 1.0
    rho_p = rho + c
    A = A0_direct - rho_p*A1
    B = A0_direct - 2*rho_p*A1 + rho_p^2*A2
    nu = A/B
    lam_hat_lg = rho_p - rho_p/(1-nu)
    lam1_lg = lam_hat_lg - c

    return Dict(
        "N"=>N, "Nprime"=>Nprime, "c"=>c, "mu1N"=>mu1, "mu2N"=>mu2, "rho"=>rho,
        "dofs_N"=>nx_s*ny_s*nz_s, "dofs_Nprime"=>Db,
        "A0"=>A0_direct, "A0_check"=>A0_check, "A0_consistency_err"=>abs(A0_direct-A0_check),
        "A2"=>A2, "B"=>B, "nu"=>nu, "cg_residual"=>resid, "cg_iters"=>ch.iters,
        "L1_LG"=>lam1_lg, "B_positive"=>(B>0), "nu_lt_1"=>(nu<1),
        "seconds"=>time()-t0
    )
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    # CLI: julia h2plus_lg_bounds.jl N Nprime c
    args = ARGS
    N = length(args)>0 ? parse(Int, args[1]) : 100
    Nprime = length(args)>1 ? parse(Int, args[2]) : N
    c = length(args)>2 ? parse(Float64, args[3]) : 1.0
    res = H2plusLG.lg_bound(N, Nprime, c)
    for (k,v) in res
        println(k, " = ", v)
    end
end
