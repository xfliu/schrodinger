# src/verified_ldl.jl
#
# Bunch-Kaufman LDLᵀ factorization of the *midpoint* of a (real or interval)
# symmetric matrix `M`, paired with a verified interval residual matrix
# bounding `M - P L D Lᵀ Pᵀ`.
#
# Replaces inline MATLAB calls in `veigs.m` (lines ~167-185, ~322-326,
# ~400-414) and `veig.m` (lines ~119-120, ~151-153, ~169-170, ~205-206).
#
# CONTRACT (soundness):
#   - When `ok == true`, every concrete matrix in the interval enclosure of
#     `ΔM` is a true residual `M_concrete - P L D Lᵀ Pᵀ` for some realisation
#     of `M_concrete ∈ M`. Equivalently, `M ⊆ P L D Lᵀ Pᵀ + ΔM` (interval
#     containment).
#   - When `ok == false`, the midpoint Bunch-Kaufman failed (e.g. singular up
#     to LAPACK tolerance). The other return values are placeholder; callers
#     should treat the result as a fall-through (D-014 spirit).
#
# Per D-005, default `uplo = :L` to match MATLAB `[L,D,P] = ldl(M)`. The `:U`
# branch exists only so the `test_verified_ldl.jl` regression can pin both
# orientations.

using LinearAlgebra
using SparseArrays
using IntervalArithmetic: Interval, interval, radius
import LDLFactorizations

"""
    verified_ldl(M; uplo=:L, ldl_backend=:auto) -> (L, D, p, ΔM, ok)

Factor the *midpoint* of `M` (a real or interval symmetric matrix) into block
LDLᵀ form via Bunch-Kaufman, and return a verified interval residual.

`ldl_backend` selects the LDL kernel:
- `:auto` (default) — sparse no-pivot LDL via LDLFactorizations.jl on
  sparse non-tridiag inputs; tridiag fast path otherwise; dense
  Bunch-Kaufman if both fail. Best for most cases.
- `:dense` — force dense Bunch-Kaufman (legacy / debugging).
- `:ldlfac` — force LDLFactorizations.jl; throws `VeigsLDLFailureError`
  on zero pivot (no dense fallback).
- `:mumps` — placeholder for MUMPS.jl integration (not wired in core;
  see seam comments in src/verified_ldl.jl). Throws on selection.

Returns:
- `L::Matrix{Float64}` — the float L factor of `mid.(M)` (lower-triangular
  when `uplo=:L`, upper-triangular when `uplo=:U`).
- `D::Matrix{Float64}` — the float block-diagonal D factor (1×1 and 2×2
  blocks; suitable input to [`inertia`](@ref)).
- `p::Vector{Int}` — the permutation vector. The float identity is
  `mid.(M)[p, p] ≈ L * D * L'` (irrespective of `uplo`; Julia's
  `BunchKaufman` returns `factor * D * factor'` either way).
- `ΔM::Matrix{Interval{Float64}}` — verified interval residual in the
  *original* ordering: `ΔM = M - P * L * D * L' * P'`, where `P` is the
  permutation matrix induced by `p`. Equivalently, `ΔM[p, p] = M[p, p] - L * D * L'`
  (interval-arithmetic).
- `ok::Bool` — `true` for any factorization LAPACK could compute, including
  the singular case (`F.info > 0` simply marks a zero pivot; the factors are
  still mathematically correct, with `D` containing a zero block). Only
  flips to `false` if LAPACK rejected the input outright (`F.info < 0`),
  which signals a programming bug rather than a numerical failure.

`uplo=:L` matches MATLAB `[L,D,P] = ldl(M)` (D-005). Sparse `M` is
densified before factoring (Julia's `bunchkaufman` is dense-only; mirrors the
MATLAB `full(mid(B))` pattern in the cluster-expansion code).

Throws `DomainError` if `mid.(M)` contains `Inf` or `NaN`.
Throws `VeigsSizeError` if `M` is non-square.

# Examples
```julia
julia> M = [4.0 2.0 1.0; 2.0 5.0 3.0; 1.0 3.0 6.0];

julia> L, D, p, ΔM, ok = verified_ldl(M);

julia> ok
true

julia> M[p, p] ≈ L * D * L'
true
```
"""
function verified_ldl(M::AbstractMatrix; uplo::Symbol = :L,
                      ldl_backend::Symbol = :auto)
    n = size(M, 1)
    n == size(M, 2) || throw(VeigsSizeError("verified_ldl: matrix must be square"))
    uplo === :L || uplo === :U ||
        throw(ArgumentError("verified_ldl: uplo must be :L or :U, got $uplo"))
    ldl_backend ∈ (:auto, :dense, :ldlfac, :mumps) ||
        throw(ArgumentError(
            "verified_ldl: ldl_backend must be :auto, :dense, :ldlfac, or :mumps (got $ldl_backend)."))
    if ldl_backend === :mumps
        throw(VeigsLDLFailureError(
            "verified_ldl: ldl_backend=:mumps is exposed but not implemented in core Veigs.jl. " *
            "MUMPS is primarily a solver — extracting L/D/p as sparse matrices for the verified " *
            "residual requires a custom MUMPS.jl integration (likely as a package extension). " *
            "See _verified_ldl_mumps_stub in src/verified_ldl.jl for the seam. " *
            "For now, use ldl_backend=:auto (defaults to LDLFactorizations.jl + dense fallback)."))
    end

    is_interval_input = M isa AbstractMatrix{<:Interval}

    # Defensive: midpoint must be finite (matches D-015 spirit)
    Mm = is_interval_input ? mid.(M) : M
    if any(!isfinite, Mm)
        throw(DomainError(M, "verified_ldl: midpoint contains Inf or NaN"))
    end

    # Tridiagonal fast path — for inputs whose midpoint has bandwidth ≤ 1
    # (1D Laplacians, banded n=2 problems, FEM 1D-line). LDL of a
    # SymTridiagonal is O(n); the full residual L·D·Lᵀ is pentadiagonal
    # so we only build O(n) of its entries.
    if uplo === :L && ldl_backend !== :dense && _is_tridiagonal_midpoint(Mm, n)
        Mi_for_tridiag = is_interval_input ? M : interval.(M)
        result = _verified_ldl_tridiag(Mi_for_tridiag, Mm, n)
        result === nothing || return result
        # ldlt(::SymTridiagonal) hit a zero pivot — fall through to dense.
    end

    # Sparse no-pivot LDL fast path (v17, audit #7). For sparse non-tridiag
    # midpoints (2D Laplacians, FEM mass matrices, banded structures)
    # LDLFactorizations.jl provides sparse symbolic LDL with AMD fill-
    # reducing ordering — typically O(N^{1.5}) work and memory on banded
    # problems vs the dense path's O(N^3) / O(N^2). It does NOT pivot on
    # the fly, so a zero pivot causes failure; we detect that and fall
    # through to dense Bunch-Kaufman.
    if uplo === :L && (ldl_backend === :auto || ldl_backend === :ldlfac) &&
       Mm isa AbstractSparseMatrix
        result = _verified_ldl_sparse_nopivot(M, Mm, n, is_interval_input)
        if result !== nothing
            return result
        end
        # Sparse LDL failed (zero pivot or LDLFactorizations exception).
        if ldl_backend === :ldlfac
            throw(VeigsLDLFailureError(
                "verified_ldl: LDLFactorizations encountered a zero pivot, but " *
                "ldl_backend=:ldlfac was forced. Use ldl_backend=:auto to fall " *
                "back to dense Bunch-Kaufman, or :dense to skip the sparse path."))
        end
    end

    # Densify the midpoint — bunchkaufman is dense-only.
    Mm_dense = Mm isa AbstractSparseMatrix ? Matrix(Mm) : Matrix(Mm)

    # Defensive symmetrization on the midpoint (cheap; same convention as
    # verified_isspd). The interval residual ΔM still bounds the full
    # asymmetry by including everything in M.
    Mm_sym = (Mm_dense + Mm_dense') / 2

    F = bunchkaufman(Symmetric(Mm_sym, uplo); check = false)

    # LAPACK signals: info < 0 → invalid argument (true failure);
    #                 info > 0 → factorization computed but a pivot is zero
    #                            (singular input — D has a zero block but the
    #                             factor is still mathematically valid).
    # MATLAB's `ldl` similarly succeeds on singular inputs, so we mirror.
    if F.info < 0
        Mi_for_failed = is_interval_input ? M : interval.(M)
        return _failed_ldl(n, Mi_for_failed)
    end

    # Extract factor as Float64 dense matrix.
    factor = uplo === :L ? Matrix(F.L) : Matrix(F.U)
    D      = Matrix(F.D)
    p      = collect(F.p)

    # Reconstruction in permuted basis: factor * D * factor' (Julia's
    # BunchKaufman convention; holds for both uplo values).
    #
    # We deliberately compute the reconstruction in *float* arithmetic
    # (BLAS gemm) and bound the rounding error with a Wilkinson-style
    # per-entry bound. The naive form `interval.(factor) * interval.(D) *
    # transpose(interval.(factor))` would do an O(n³) interval matmul with
    # no BLAS path — ~100× slower and the bottleneck for n ≥ 100.
    #
    # Backward-error bound (Higham, "Accuracy & Stability", §10):
    #   |fl(L D Lᵀ) − exact(L D Lᵀ)|_ij  ≤  γ_{3n} · (|L|·|D|·|Lᵀ|)_ij
    # with γ_k = k·u / (1 − k·u), u = unit roundoff = eps()/2.
    recon_f = factor * D * transpose(factor)        # BLAS, O(n³)

    abs_L   = abs.(factor)
    abs_D   = abs.(D)
    LD_abs  = abs_L * abs_D                         # |L|·|D|
    u       = eps(Float64) / 2
    k       = 3 * n
    γ       = k * u >= 1.0 ? Inf : k * u / (1 - k * u)
    B_round = LD_abs * transpose(abs_L)             # |L|·|D|·|Lᵀ|, entry-wise
    B_round .= γ .* B_round                         # γ_{3n}·(|L|·|D|·|Lᵀ|)_ij

    # ΔM[p, p] = M[p, p] − exact(L D Lᵀ). Compute as a (mid, rad) split in
    # float, then pack to `Matrix{Interval{Float64}}` once at the end via
    # outward-rounded addition (`_pack_midrad` from src/lehmann_behnke.jl).
    # This replaces O(n²) per-entry interval-arithmetic ops with two float
    # matrix subtractions and one packing pass — BLAS-friendly and
    # ~10× faster than the per-entry interval form at n ≥ 1000.
    #
    # Soundness: ΔM_mid = fl(Mm[p,p] − recon_f) introduces rounding
    # bounded by u·max(|Mm[p,p]|, |recon_f|). Combined with input radius
    # (if any) and the Wilkinson bound, the packed interval rigorously
    # encloses M[p,p] − exact(L·D·Lᵀ).
    Mm_perm  = Mm_dense[p, p]
    ΔM_mid   = Mm_perm .- recon_f                   # float sub, dense
    sub_slack = u .* (abs.(Mm_perm) .+ abs.(recon_f))   # float-sub rounding bound
    ΔM_rad   = B_round .+ sub_slack
    if is_interval_input
        # Add input radius (densify if sparse — same order as Mm_dense).
        Mr_perm = if M isa AbstractSparseMatrix
            Matrix(radius.(M))[p, p]
        else
            radius.(M)[p, p]
        end
        ΔM_rad .+= Mr_perm
    end

    ΔM_perm = _pack_midrad(ΔM_mid, ΔM_rad, Float64)

    # Un-permute to original index space.
    invp = invperm(p)
    ΔM   = ΔM_perm[invp, invp]

    return (factor, D, p, ΔM, true)
end

# Placeholder return when the midpoint factorization fails.
function _failed_ldl(n::Int, Mi::AbstractMatrix)
    L  = zeros(n, n)
    D  = zeros(n, n)
    p  = collect(1:n)
    ΔM = Matrix(Mi)            # caller treats as garbage when ok == false
    return (L, D, p, ΔM, false)
end

# ============================================================================
# Tridiagonal fast path
# ----------------------------------------------------------------------------
# For an n×n symmetric tridiagonal A (bandwidth 1), the LDLᵀ factorization
# preserves bandwidth: L is unit lower bidiagonal and D is purely diagonal
# (no 2×2 blocks). Both factor and reconstruction are O(n) work.
#
# `LinearAlgebra.ldlt(::SymTridiagonal)` does the unpivoted LDL in stdlib;
# it can fail with `ZeroPivotException` on rank-deficient indefinite
# inputs (e.g., when M = A − λB and λ sits exactly on an eigenvalue).
# We catch that case and signal the caller to fall through to the dense
# Bunch-Kaufman path.
#
# This brings n=5000 sparse-Laplacian `verified_ldl` from ~12 s (dense BK
# + 3 dense gemms) to ~ms.
# ============================================================================

function _is_tridiagonal_midpoint(Mm::AbstractMatrix, n::Int)
    n ≥ 2 || return true                    # 1×1 and 0×0 trivially tridiag
    @inbounds for j in 1:n
        for i in 1:n
            if abs(i - j) > 1 && !iszero(Mm[i, j])
                return false
            end
        end
    end
    return true
end

# Fast O(nnz) variant for sparse matrices. Walks the CSC nzval directly via
# colptr/rowval — no scalar `getindex` (which would be O(log nnz) per probe
# and turns the dense version into O(n² log nnz) overall).
function _is_tridiagonal_midpoint(Mm::SparseMatrixCSC, n::Int)
    n ≥ 2 || return true
    rowval = SparseArrays.rowvals(Mm)
    colptr = SparseArrays.getcolptr(Mm)
    @inbounds for j in 1:n
        for k in colptr[j]:(colptr[j+1] - 1)
            i = rowval[k]
            if abs(i - j) > 1 && !iszero(Mm.nzval[k])
                return false
            end
        end
    end
    return true
end

function _verified_ldl_tridiag(Mi::AbstractMatrix, Mm::AbstractMatrix, n::Int)
    # Extract main diagonal and (symmetrized) sub-diagonal.
    dv = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        dv[i] = Float64(Mm[i, i])
    end
    ev = Vector{Float64}(undef, max(n - 1, 0))
    @inbounds for i in 1:(n - 1)
        ev[i] = 0.5 * (Float64(Mm[i + 1, i]) + Float64(Mm[i, i + 1]))
    end

    T = SymTridiagonal(dv, ev)
    F = try
        ldlt(T)
    catch err
        err isa LinearAlgebra.ZeroPivotException || rethrow(err)
        return nothing
    end

    # `ldlt(::SymTridiagonal)` returns LDLt{T,SymTridiagonal{T,V}} whose
    # `.data.dv` is D's diagonal and whose `.data.ev` is L's subdiagonal
    # (with implicit unit diagonal on L).
    D_diag = F.data.dv
    L_sub  = F.data.ev

    if any(!isfinite, D_diag) || any(!isfinite, L_sub)
        return nothing
    end

    # Return L as a sparse bidiagonal matrix and D as a `Diagonal{Float64}`
    # rather than dense n×n. This keeps memory at O(n) and lets `inertia(D)`
    # skip the O(n²) `interval.(D)` broadcast over the all-zero off-diagonal.
    # `inertia` already has a `::Diagonal` fast path; downstream callers
    # only use `D` via `inertia` and `L` via norm queries on the residual,
    # so the type narrowing is transparent.
    L = sparse(SparseArrays.spdiagm(0 => ones(Float64, n),
                                    -1 => copy(L_sub)))
    D = Diagonal(copy(D_diag))
    p = collect(1:n)        # SymTridiagonal LDL has no pivoting

    # Build the residual `M − L·D·Lᵀ` directly. L is unit bidiag so
    # `L·D·Lᵀ` is pentadiagonal — only entries within bandwidth 2 of the
    # diagonal are nonzero. We compute those entries closed-form to avoid
    # an O(n²) gemm.
    #
    # With L[i,i] = 1 and L[i+1,i] = e_i (subdiag), and D diagonal:
    #     (L·D·Lᵀ)[i,j] = Σ_k L[i,k] · D[k,k] · L[j,k]
    # The only nonzero L[i,k] are k = i (= 1) and k = i − 1 (= e_{i-1}).
    #
    # Wilkinson per-entry bound (Higham §10):
    #     |fl(L D Lᵀ) − exact|_{ij} ≤ γ_{3n} · (|L| |D| |Lᵀ|)_{ij}
    # which is the same closed form with absolute values.
    u   = eps(Float64) / 2
    k_w = 3 * n
    γ   = k_w * u >= 1.0 ? Inf : k_w * u / (1 - k_w * u)

    function L_at(i, k)
        i == k && return 1.0
        i == k + 1 && k ≥ 1 && k ≤ n - 1 && return L_sub[k]
        return 0.0
    end

    # ΔM = Mi − exact(L D Lᵀ). The reconstruction L·D·Lᵀ is pentadiagonal
    # (bandwidth 2), and Mi is tridiagonal (bandwidth 1) — so ΔM is also
    # pentadiagonal at most. Build it as a `SparseMatrixCSC{Interval{Float64}}`
    # with at most 5n nonzeros instead of densifying to ~n² entries
    # (~400 MB at n=5000) — DESIGN.md §1.
    nz_estimate = min(5 * n, n * n)
    rows  = Vector{Int}(undef, 0); sizehint!(rows, nz_estimate)
    cols  = Vector{Int}(undef, 0); sizehint!(cols, nz_estimate)
    vals  = Vector{Interval{Float64}}(undef, 0); sizehint!(vals, nz_estimate)

    @inbounds for j in 1:n
        for i in max(1, j - 2):min(n, j + 2)
            recon = 0.0
            babs  = 0.0
            for κ in max(1, max(i, j) - 1):min(n, min(i, j))
                Lik = L_at(i, κ)
                Ljk = L_at(j, κ)
                Lik == 0.0 && continue
                Ljk == 0.0 && continue
                recon += Lik * D_diag[κ] * Ljk
                babs  += abs(Lik) * abs(D_diag[κ]) * abs(Ljk)
            end
            b = γ * babs
            slack = isfinite(b) ? interval(-b, b) : interval(-Inf, Inf)
            entry = Mi[i, j] - interval(recon) + slack
            push!(rows, i); push!(cols, j); push!(vals, entry)
        end
    end

    ΔM = sparse(rows, cols, vals, n, n)
    return (L, D, p, ΔM, true)
end

# ============================================================================
# Sparse no-pivot LDL fast path (v17, audit #7)
# ----------------------------------------------------------------------------
# `LDLFactorizations.jl` provides a pure-Julia sparse symbolic-then-numerical
# LDLᵀ factorization with AMD fill-reducing ordering. It does NOT pivot on
# the fly: a zero pivot causes failure (returned as a zero in `F.d`). For
# generic indefinite matrices we detect that case and signal the caller to
# fall through to dense Bunch-Kaufman.
#
# When successful, this brings 2D Laplacian P1 FEM verified_ldl from O(N³)
# dense Bunch-Kaufman down to O(N^{1.5}) sparse LDL — orders of magnitude
# faster at N ≥ 1000 (DESIGN.md §1).
#
# Per D-013, the residual matrix returned is a SparseMatrixCSC{Interval}.
# ============================================================================

function _verified_ldl_sparse_nopivot(M::AbstractMatrix, Mm::SparseMatrixCSC,
                                      n::Int, is_interval_input::Bool)
    # Symmetrize the sparse midpoint (cheap, same convention as dense path).
    Mm_sym = (Mm + transpose(Mm)) / 2

    # Try sparse LDL. LDLFactorizations may either throw or signal failure
    # via zero entries in `F.d` — handle both.
    F = try
        LDLFactorizations.ldl(Mm_sym)
    catch
        return nothing
    end

    d = F.d
    L_sparse = F.L                  # sparse unit lower triangular, no diag stored
    p        = collect(F.P)         # AMD permutation

    # Detect zero / non-finite pivots (factorization failed or matrix is
    # singular). LDLFactorizations leaves `d[i] = 0` (or garbage) past a
    # zero pivot — fall through to dense rather than risk garbage outputs.
    @inbounds for i in 1:n
        di = d[i]
        if !isfinite(di) || di == 0.0
            return nothing
        end
    end

    # `F.L` is unit lower triangular without the diagonal stored. Add
    # the diagonal so downstream consumers see a complete L (matches the
    # dense path's `Matrix(F.L)` convention, including the unit diagonal).
    L = L_sparse + sparse(I, n, n)
    D = Diagonal(copy(d))

    # Float reconstruction `L · D · Lᵀ` (sparse matmul, O(nnz(L)·nnz(L)/n)).
    DLt    = D * transpose(L)
    recon_f = L * DLt                                        # sparse N×N

    # Wilkinson per-entry bound on the reconstruction (sparse matmul too).
    abs_L  = abs.(L)
    abs_D  = Diagonal(abs.(d))
    LD_abs = abs_L * abs_D
    B_round = LD_abs * transpose(abs_L)
    u  = eps(Float64) / 2
    kw = 3 * n
    γ  = kw * u >= 1.0 ? Inf : kw * u / (1 - kw * u)
    B_round = γ .* B_round

    # ΔM[p, p] = M[p, p] − fl(L · D · Lᵀ)
    Mm_perm = Mm[p, p]
    ΔM_mid  = Mm_perm - recon_f

    # Float-subtraction rounding bound (entry-wise): u·(|a| + |b|).
    sub_slack = u .* (abs.(Mm_perm) .+ abs.(recon_f))
    ΔM_rad    = B_round + sub_slack

    if is_interval_input
        # Add the input's per-entry radius (sparse if M sparse).
        Mr_perm = if M isa SparseMatrixCSC
            M_radius_sparse = SparseMatrixCSC(M.m, M.n,
                                              copy(SparseArrays.getcolptr(M)),
                                              copy(SparseArrays.rowvals(M)),
                                              Float64.(radius.(nonzeros(M))))
            M_radius_sparse[p, p]
        else
            sparse(radius.(M))[p, p]
        end
        ΔM_rad = ΔM_rad + Mr_perm
    end

    # Pack as SparseMatrixCSC{Interval{Float64}} via outward-rounded
    # addition. The result's pattern is the union of (Mm_perm, recon_f,
    # B_round) — typically dominated by `B_round`'s symbolic pattern
    # (|L|·|D|·|Lᵀ|), which has the symbolic fill of L.
    ΔM_perm = _pack_midrad_sparse(ΔM_mid, ΔM_rad)

    # Un-permute to original index space.
    invp = invperm(p)
    ΔM   = ΔM_perm[invp, invp]

    return (L, D, p, ΔM, true)
end

# Sparse counterpart to `_pack_midrad` (defined in src/lehmann_behnke.jl).
# Builds a SparseMatrixCSC{Interval{Float64}} where each (i,j) entry is
# `interval(mid_mat[i,j]) + interval(-rad_mat[i,j], rad_mat[i,j])` with
# directed-rounded addition. The sparsity pattern is the union of the
# nonzero patterns of mid and rad (so a position with mid=0 but rad>0
# becomes an `interval(-rad, rad)` entry).
# ============================================================================
# Seam for option 2: MUMPS.jl integration (not implemented in core)
# ----------------------------------------------------------------------------
# MUMPS provides high-performance sparse symmetric indefinite LDL with
# 1×1 / 2×2 pivoting (the equivalent of MA57 used by MATLAB+INTLAB). It
# would close the remaining 1.3–1.5× gap to MATLAB on 2D problems, and
# gracefully handle indefinite cases where LDLFactorizations.jl hits a
# zero pivot (currently we fall back to dense Bunch-Kaufman in that case).
#
# Why it isn't wired by default:
#   1. MUMPS is primarily a SOLVER — it computes Ax = b in internal
#      state, and exposing L, D, p as `SparseMatrixCSC` requires
#      navigating MUMPS_jll's Fortran-layout output. The current
#      verified_ldl interface returns (L, D, p, ΔM, ok); a clean MUMPS
#      backend likely calls for a refactor toward "inertia(M − λ B) at
#      λ" as the primitive (which MUMPS exposes via INFOG(12)).
#   2. MUMPS.jl is a heavy dep (C++ MUMPS_jll binary, fragile build on
#      some macOS/Apple Silicon configurations).
#
# How to wire it (suggested approach via Julia package extensions):
#
#   1. Add to Project.toml:
#        [weakdeps]
#        MUMPS = "..."
#        [extensions]
#        VeigsMUMPSExt = "MUMPS"
#
#   2. Create ext/VeigsMUMPSExt.jl with:
#        module VeigsMUMPSExt
#        using Veigs, MUMPS, SparseArrays, IntervalArithmetic
#        function Veigs._verified_ldl_mumps_impl(M, Mm, n, is_interval_input)
#            # ... use MUMPS.jl to factorize, extract L/D/p, build ΔM
#        end
#        end
#
#   3. Replace the throw at the top of verified_ldl with a dispatch:
#        if ldl_backend === :mumps
#            return _verified_ldl_mumps_impl(M, Mm, n, is_interval_input)
#        end
#      The default (no extension loaded) implementation throws the
#      VeigsLDLFailureError currently in place.
#
# Until this is implemented, ldl_backend=:mumps throws an informative
# error directing users here. Contributions welcome.
# ============================================================================

function _pack_midrad_sparse(mid_mat::SparseMatrixCSC{<:Real},
                             rad_mat::SparseMatrixCSC{<:Real})
    n_rows, n_cols = size(mid_mat)
    size(rad_mat) == (n_rows, n_cols) ||
        throw(DimensionMismatch("_pack_midrad_sparse: mid and rad shape mismatch"))

    # Union pattern: build a 1-valued sparse mask combining both patterns.
    pattern = (mid_mat .!= 0) .| (rad_mat .!= 0)
    rows_p, cols_p, _ = findnz(pattern)

    rows = Vector{Int}(undef, length(rows_p))
    cols = Vector{Int}(undef, length(cols_p))
    vals = Vector{Interval{Float64}}(undef, length(rows_p))
    @inbounds for k in eachindex(rows_p)
        i = rows_p[k]; j = cols_p[k]
        m = Float64(mid_mat[i, j])
        r = Float64(rad_mat[i, j])
        rows[k] = i; cols[k] = j
        vals[k] = if isfinite(r)
            interval(m) + interval(-r, r)
        else
            interval(-Inf, Inf)
        end
    end
    return sparse(rows, cols, vals, n_rows, n_cols)
end
