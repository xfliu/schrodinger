# src/lehmann_behnke.jl
#
# Port of MATLAB `Lehmann_Behnke__` (subfunction in `veigs.m`, lines ~575-635).
# Verified bounds for an eigenvalue cluster (λ_r, …, λ_s) of the generalized
# symmetric eigenproblem `A x = λ B x`, using complementary variational
# principles with separators ρ (upper bound on λ_{r-1}) and σ (lower bound
# on λ_{s+1}).
#
# CONTRACT (soundness):
#   The returned `Vector{Interval{Float64}}` of length `s - r + 1` contains
#   verified enclosures for the cluster eigenvalues, after un-doing the
#   optional shift. If the Lehmann SB matrix fails the SPD test on either
#   side, `VeigsLehmannBehnkeError` is thrown rather than returning a wrong
#   answer.
#
# Algorithm sketch:
#   1. Optional shift `A := A − λ_shift * B`, where `λ_shift = eig_list[r]`.
#      Adjusts ρ and σ correspondingly.
#   2. A0 = V'BV, A1 = V'AV (interval m×m).
#      nV  = mid(B) \ (mid(A) * mid(V))      (float preconditioner)
#      Err = (B*nV − A*V)' (B*nV − A*V) / λ_B_min
#      A2  = V' A nV − nV' (B nV − A V) + Err.
#   3. Lower side: SA = -A2 + σ A1, SB = -A1 + σ A0. veig(SA, SB) gives `low`.
#   4. Upper side: SA = A2 - ρ A1, SB = A1 - ρ A0. veig gives `upper`.
#   5. λ = hull(low + λ_shift, upper + λ_shift)  elementwise.
#
# `eig_list` is required (not just `r, s`): the shift uses `eig_list[r]`.

using LinearAlgebra
using IntervalArithmetic: Interval, interval, mid, hull, inf, sup

"""
    lehmann_behnke(A, B, eig_list, V, ρ, σ, λ_B_min, r, s; do_shift=true)
        -> Vector{Interval{Float64}}

Verified bounds for the eigenvalue cluster `(λ_r, …, λ_s)` of `A x = λ B x`.

Arguments:
- `A`, `B::AbstractMatrix` — symmetric `n×n` (real or interval). `B` PD.
- `eig_list::AbstractVector{<:Real}` — the float-approximate eigenvalues
  (sorted ascending). `eig_list[r]` is used as the shift target when
  `do_shift=true`.
- `V::AbstractMatrix` — float `n×m` matrix of approximate eigenvectors for
  indices `r:s` (so `m = s - r + 1`). Need not be `B`-orthonormal.
- `ρ::Real` — verified upper bound on `λ_{r-1}` (left separator).
- `σ::Real` — verified lower bound on `λ_{s+1}` (right separator).
- `λ_B_min::Real` — verified lower bound on the smallest eigenvalue of `B`
  (positive). Used as the energy denominator in the residual term.
- `r, s::Int` — cluster index range (`1 ≤ r ≤ s ≤ n`).

Keyword:
- `do_shift::Bool = true` — apply the conditioning shift `A → A − eig_list[r] B`.

Returns: `Vector{Interval{Float64}}` of length `m = s - r + 1`. For `r == s`
the vector has one element (the scalar Rayleigh-style ratio).

Throws `VeigsLehmannBehnkeError` if either `SB` matrix fails to verify PD.
"""
function lehmann_behnke(A::AbstractMatrix, B::AbstractMatrix,
                        eig_list::AbstractVector{<:Real},
                        V::AbstractMatrix,
                        ρ::Real, σ::Real, λ_B_min::Real,
                        r::Int, s::Int;
                        do_shift::Bool = true,
                        precision::Type{<:AbstractFloat} = Float64,
                        outer_matmul::Symbol = :fast)
    outer_matmul ∈ (:fast, :tight) ||
        throw(ArgumentError("lehmann_behnke: outer_matmul must be :fast or :tight, got $outer_matmul"))
    n = size(A, 1)
    n == size(A, 2) || throw(VeigsSizeError("lehmann_behnke: A must be square"))
    size(B) == size(A) || throw(VeigsSizeError("lehmann_behnke: A and B must have the same shape"))
    size(V, 1) == n || throw(VeigsSizeError("lehmann_behnke: V must have $n rows, got $(size(V, 1))"))
    1 ≤ r ≤ s ≤ n || throw(ArgumentError("lehmann_behnke: need 1 ≤ r ≤ s ≤ n, got r=$r s=$s n=$n"))
    size(V, 2) == s - r + 1 ||
        throw(ArgumentError("lehmann_behnke: V must have s-r+1 = $(s-r+1) columns, got $(size(V, 2))"))
    λ_B_min > 0 ||
        throw(ArgumentError("lehmann_behnke: λ_B_min must be positive, got $λ_B_min"))

    # Promote to interval at the requested precision. For `precision === Float64`
    # this is the existing behaviour; for higher precision (`BigFloat`,
    # `Double64`, …) every subsequent op runs at that precision — including the
    # n×m outer matmuls, which is the only way the inner-step bound width
    # actually tightens (DESIGN.md §3). At n=5000 BigFloat is ~100× slower
    # than Float64; the option is intended for *demonstrating* the bound floor
    # is set by the algorithm, not the arithmetic, on small problems.
    Ai = _lb_promote_interval(A, precision)
    Bi = _lb_promote_interval(B, precision)
    Vf = V isa AbstractMatrix{<:Interval} ? precision.(mid.(V)) : convert(Matrix{precision}, V)

    # ---- Step 1: optional shift -------------------------------------------
    λ_shift = zero(precision)
    if do_shift
        λ_shift = precision(eig_list[r])
        Ai = Ai .- interval(λ_shift) .* Bi
        # Adjust separators outward (sup/inf to be conservative).
        ρ = sup(interval(precision(ρ)) - interval(λ_shift))
        σ = inf(interval(precision(σ)) - interval(λ_shift))
    end

    # ---- Step 2: Lehmann submatrices --------------------------------------
    # The four n×m outer matmuls (Bi·V, Ai·V, Bi·nV, Ai·nV) are the
    # dominant cost in this routine. Two paths, controlled by `outer_matmul`:
    #
    # `:fast` (default) — float gemm + per-entry Wilkinson scalar bound.
    #   `IntervalArithmetic.jl`'s interval×float kernel is ~1000× slower
    #   than pure float sparse·dense at n=5000. We bypass it: midpoint
    #   via float gemm, radius via |fl(M·X) − M·X|_ij ≤ γ_k · (|M|·|X|)_ij
    #   with γ_k = k·u/(1−k·u), k = max-nnz-per-row(M) (zero terms
    #   contribute exactly to the dot product, so the chain length is
    #   the actual nnz count, not the inner dim n).
    #   Trade-off: 5–10× faster than `:tight` at n=5000, but the
    #   per-entry radius is the worst-case Higham bound, ~10× looser
    #   than the actual rounding for matrices with integer multipliers
    #   (1D Laplacian, FEM mass matrices) where mults are exact.
    #
    # `:tight` — IntervalArithmetic.jl's directed-rounded interval×float
    #   matmul. Slower but per-entry intervals are tightest possible
    #   (often exact for integer-multiplier sparse problems).
    #
    # Both paths are bound-equivalent on dense problems within ULP;
    # the difference matters mainly for tridiagonal-and-similar
    # problems where individual products happen to be exact floats.
    A_mid_f = precision === Float64 ? _mid_float(Ai) : _mid_promoted(Ai, precision)
    B_mid_f = precision === Float64 ? _mid_float(Bi) : _mid_promoted(Bi, precision)

    # nV = mid(B) \ (mid(A) * mid(V))    (float preconditioner; identical in both paths)
    nV    = B_mid_f \ (A_mid_f * Vf)

    if outer_matmul === :tight
        # v12-equivalent path: directed-rounded interval × float matmul.
        # Each per-entry interval has width = actual rounding of the
        # specific dot product, often = 0 (exact) for integer-multiplier M.
        BV  = Bi * Vf
        AV  = Ai * Vf
        BnV = Bi * nV
        AnV = Ai * nV
    else
        # v14 path: float gemm + Wilkinson scalar bound.
        A_rad_f = _rad_float(Ai, precision)
        B_rad_f = _rad_float(Bi, precision)
        abs_A   = abs.(A_mid_f)
        abs_B   = abs.(B_mid_f)
        abs_Vf  = abs.(Vf)
        abs_nV  = abs.(nV)
        u       = eps(precision) / 2
        k_A = _max_nnz_per_row(A_mid_f)
        k_B = _max_nnz_per_row(B_mid_f)
        γA  = (k_A * u) / (1 - k_A * u)
        γB  = (k_B * u) / (1 - k_B * u)

        BV_mid  = B_mid_f * Vf
        AV_mid  = A_mid_f * Vf
        BnV_mid = B_mid_f * nV
        AnV_mid = A_mid_f * nV

        BV_rad  = γB .* (abs_B * abs_Vf)
        AV_rad  = γA .* (abs_A * abs_Vf)
        BnV_rad = γB .* (abs_B * abs_nV)
        AnV_rad = γA .* (abs_A * abs_nV)
        if B_rad_f !== nothing
            BV_rad  = BV_rad  .+ B_rad_f * abs_Vf
            BnV_rad = BnV_rad .+ B_rad_f * abs_nV
        end
        if A_rad_f !== nothing
            AV_rad  = AV_rad  .+ A_rad_f * abs_Vf
            AnV_rad = AnV_rad .+ A_rad_f * abs_nV
        end

        BV  = _pack_midrad(BV_mid,  BV_rad,  precision)
        AV  = _pack_midrad(AV_mid,  AV_rad,  precision)
        BnV = _pack_midrad(BnV_mid, BnV_rad, precision)
        AnV = _pack_midrad(AnV_mid, AnV_rad, precision)
    end

    # m×m interval matmuls: cheap because m ≤ 30.
    A0 = transpose(Vf) * BV                # m×m interval
    A1 = transpose(Vf) * AV                # m×m interval

    # Err = (B·nV − A·V)' (B·nV − A·V) / λ_B_min
    err_vec = BnV .- AV                    # n×m interval
    Err     = (transpose(err_vec) * err_vec) ./ interval(precision(λ_B_min))

    # A2 = V' · A · nV − nV' · (B·nV − A·V) + Err
    A2 = transpose(Vf) * AnV .-
         transpose(nV) * err_vec .+
         Err

    # ---- Step 3: lower bound (σ side) -------------------------------------
    σ_int = interval(precision(σ))
    SA_lo = .-A2 .+ σ_int .* A1
    SB_lo = .-A1 .+ σ_int .* A0

    verified_isspd(sym_hull(SB_lo)) ||
        throw(VeigsLehmannBehnkeError(
            "lehmann_behnke: SB (lower-side) is not verified PD; cluster too wide or σ too tight"))

    low = _lehmann_solve(SA_lo, SB_lo, r, s)

    # ---- Step 4: upper bound (ρ side) -------------------------------------
    ρ_int = interval(precision(ρ))
    SA_up = A2 .- ρ_int .* A1
    SB_up = A1 .- ρ_int .* A0

    verified_isspd(sym_hull(SB_up)) ||
        throw(VeigsLehmannBehnkeError(
            "lehmann_behnke: SB (upper-side) is not verified PD; cluster too wide or ρ too tight"))

    upper = _lehmann_solve(SA_up, SB_up, r, s)

    # ---- Step 5: compose verified interval (re-shift if needed) -----------
    @assert length(low) == length(upper) == s - r + 1
    λ_shift_int = interval(λ_shift)
    out = Vector{Interval{precision}}(undef, s - r + 1)
    for i in eachindex(out)
        l = low[i] + λ_shift_int
        u = upper[i] + λ_shift_int
        out[i] = interval(inf(l), sup(u))
    end
    return out
end

# Promote a (real or interval) matrix to `Matrix{Interval{T}}` (or
# `SparseMatrixCSC{Interval{T}, …}` when sparse), preserving sparsity.
# Returns the input unchanged when its element type is already
# `Interval{T}`. Used by `lehmann_behnke` precision option (DESIGN.md §3).
function _lb_promote_interval(M::AbstractMatrix, ::Type{T}) where T<:AbstractFloat
    if M isa AbstractMatrix{<:Interval}
        if eltype(eltype(M)) === T
            return M
        end
        # Element-wise widen / narrow each interval.
        return _convert_interval_eltype.(M, T)
    end
    if M isa AbstractSparseMatrix
        return interval.(T.(M))
    end
    return interval.(convert(Matrix{T}, M))
end

@inline _convert_interval_eltype(x::Interval, ::Type{T}) where T<:AbstractFloat =
    interval(T(inf(x)), T(sup(x)))

# Extract the float midpoint (eltype Float64). Sparsity-preserving on
# SparseMatrixCSC. Used by the v14 outer-matmul refactor.
_mid_float(M::AbstractMatrix{<:Interval}) = mid.(M)
_mid_float(M::AbstractMatrix) = M

# Extract the midpoint at a higher precision (`T <: AbstractFloat`).
function _mid_promoted(M::AbstractMatrix, ::Type{T}) where T<:AbstractFloat
    if M isa AbstractMatrix{<:Interval}
        return T.(mid.(M))
    end
    return T.(M)
end

# Extract the radius matrix as float, or `nothing` if M is float-only
# (radius is identically zero — saves a matmul in the radius bound).
function _rad_float(M::AbstractMatrix, ::Type{T}) where T<:AbstractFloat
    if M isa AbstractMatrix{<:Interval}
        return T.(radius.(M))
    end
    return nothing
end

# Maximum number of nonzero entries in any row of `M`. For sparse M the
# Wilkinson chain length on `(M·X)[i,j]` equals this count (zero terms
# contribute exactly zero). For dense M the chain length is the full
# inner dimension.
_max_nnz_per_row(M::AbstractMatrix) = size(M, 2)
function _max_nnz_per_row(M::SparseMatrixCSC)
    n = size(M, 1)
    n == 0 && return 0
    counts  = zeros(Int, n)
    rowval  = SparseArrays.rowvals(M)
    @inbounds for k in eachindex(rowval)
        counts[rowval[k]] += 1
    end
    return maximum(counts)
end

# Build `Matrix{Interval{T}}` from float (mid, rad) pair via outward-
# rounded interval arithmetic. The `interval(m) + interval(-r, r)` form
# delegates to IntervalArithmetic.jl's directed-rounding addition, so
# the resulting endpoints rigorously enclose the true mid ± rad value
# even if the float ± would have rounded inward.
function _pack_midrad(mid_mat::AbstractMatrix, rad_mat::AbstractMatrix,
                      ::Type{T}) where T<:AbstractFloat
    out = Matrix{Interval{T}}(undef, size(mid_mat))
    @inbounds for j in 1:size(mid_mat, 2)
        for i in 1:size(mid_mat, 1)
            r = rad_mat[i, j]
            out[i, j] = interval(T(mid_mat[i, j])) + interval(-T(r), T(r))
        end
    end
    return out
end

# Solve the small Lehmann subproblem and return a vector of m intervals.
# Scalar case (s == r) is a 1×1 ratio; cluster case calls `veig`.
function _lehmann_solve(SA::AbstractMatrix{<:Interval},
                        SB::AbstractMatrix{<:Interval}, r::Int, s::Int)
    m = s - r + 1
    if m == 1
        # Scalar ratio: SA[1,1] / SB[1,1]
        return [SA[1, 1] / SB[1, 1]]
    end
    # Cluster case: use veig on the symmetrized small EVP.
    bounds, _ind = veig(sym_hull(SA), sym_hull(SB))
    # `veig` may return more rows than requested if eigenvalues cluster, but
    # we always need exactly `m` outputs. Take the first `m` (matches MATLAB
    # default output of `veig` on full range).
    length(bounds) ≥ m ||
        throw(VeigsLehmannBehnkeError(
            "lehmann_behnke: veig returned $(length(bounds)) bounds, expected ≥ $m"))
    return bounds[1:m]
end
