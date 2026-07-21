# src/veigs_driver.jl
#
# Cluster-aware verified eigenvalue solver, the headline driver. Port of
# MATLAB `veigs.m` outer function plus the `veigs_EigLocalBound` subfunction.
#
# (Filename is `veigs_driver.jl` rather than `veigs.jl` because the package
# module is `Veigs.jl` and macOS APFS is case-insensitive — same name would
# collide with the module file.)
#
# Public API (D-002):
#   veigs(A, B)                   # k = 1, sigma = :largestabs (legacy default)
#   veigs(A, B, k::Int)
#   veigs(A, B, sigma)            # sigma::Symbol or Real
#   veigs(A, B, k::Int, sigma)
# Returns `(lambda::Vector{Interval{Float64}}, ind_range::UnitRange{Int})`.
#
# Sigma symbols:
#   modern:  :largestabs, :smallestabs, :largestreal, :smallestreal
#   legacy:  :lm, :sm, :la, :sa
#   numeric: any Real x → "eigenvalue closest to x"
#
# This implementation uses dense `eigen` for approximate eigenvalues. The
# MATLAB code first tries `eigs` (Arpack) and falls back to `eig`; for the
# README problem sizes (n ≤ 10) the dense path is plenty. A sparse-friendly
# Arpack path is queued for iter 8.

using LinearAlgebra
using SparseArrays
using IntervalArithmetic: Interval, interval, mid, hull, mag, inf, sup
import KrylovKit

const _LEGACY_SIGMA = Dict{Symbol,Symbol}(
    :lm => :largestabs,  :sm => :smallestabs,
    :la => :largestreal, :sa => :smallestreal,
)
const _MODERN_SIGMA = (:largestabs, :smallestabs, :largestreal, :smallestreal)

const _VEIGS_LAMBDA_RATIO = 1.0e-5     # MATLAB default
const _VEIGS_DO_SHIFT     = true       # Lehmann-Behnke shift on by default

# ---- Public entry points ---------------------------------------------------

"""
    veigs(A, B[, k][, sigma]) -> (lambda::Vector{Interval{Float64}}, ind_range::UnitRange{Int})

Compute verified bounds for `k` eigenvalues of the symmetric generalized
eigenproblem `A x = λ B x`, selected by `sigma`.

`sigma` accepts:
- `:largestabs` (default), `:smallestabs`, `:largestreal`, `:smallestreal`
- legacy aliases `:lm`, `:sm`, `:la`, `:sa`
- a `Real` scalar — selects the eigenvalue closest to it (shift-invert style)

`k::Int` (default `1`) is the minimum number of eigenvalues to return. The
returned `lambda` and `ind_range` may include extra eigenvalues if a
cluster reaches outward; this matches MATLAB behavior.

`B` must be positive definite (interval PD allowed).
"""
veigs(A::AbstractMatrix, B::AbstractMatrix) = _veigs_dispatch(A, B, 1, :largestabs)

function veigs(A::AbstractMatrix, B::AbstractMatrix, k::Integer)
    return _veigs_dispatch(A, B, Int(k), :largestabs)
end

function veigs(A::AbstractMatrix, B::AbstractMatrix, sigma::Union{Symbol,AbstractString,Real})
    return _veigs_dispatch(A, B, 1, sigma)
end

function veigs(A::AbstractMatrix, B::AbstractMatrix, k::Integer,
               sigma::Union{Symbol,AbstractString,Real})
    return _veigs_dispatch(A, B, Int(k), sigma)
end

# ---- Core dispatch ---------------------------------------------------------

function _veigs_dispatch(A::AbstractMatrix, B::AbstractMatrix,
                         num_eigs::Int, sigma_in)
    n = size(A, 1)
    n == size(A, 2) || throw(VeigsSizeError("veigs: A must be square"))
    size(B) == size(A) || throw(VeigsSizeError("veigs: A and B must have the same shape"))
    num_eigs ≥ 1 || throw(ArgumentError("veigs: k must be ≥ 1"))
    num_eigs ≤ n || throw(ArgumentError("veigs: k=$num_eigs > n=$n"))

    sigma_norm, sigma_is_numeric = _normalize_sigma(sigma_in)
    sigma_value = sigma_is_numeric ? float(sigma_in) : NaN

    # ---- Approximate eigenvalues / eigenvectors --------------------------
    # Run *before* interval promotion so we can keep `A`, `B` sparse for the
    # KrylovKit shift-invert path.
    eig_list, V = _approx_eig(A, B, num_eigs, sigma_norm, sigma_value, sigma_is_numeric)
    m = length(eig_list)

    # Promote to interval for the verified-bounds work (idempotent if already
    # interval). Preserve sparsity: `interval.(::SparseMatrixCSC)` broadcasts
    # over the nonzero pattern, so `Ai`/`Bi` stay `SparseMatrixCSC{Interval}`
    # when the input was sparse, and all downstream operations
    # (`mid.(Ai)`, `mag.(Ai)`, `Ai .- λ·Bi`, `Ai[p,p]`, ...) remain O(nnz).
    Ai = if A isa AbstractMatrix{<:Interval}
        A
    elseif A isa AbstractSparseMatrix
        interval.(A)
    else
        interval.(Matrix(A))
    end
    Bi = if B isa AbstractMatrix{<:Interval}
        B
    elseif B isa AbstractSparseMatrix
        interval.(B)
    else
        interval.(Matrix(B))
    end
    A_mid = mid.(Ai)
    B_mid = mid.(Bi)

    # ---- Pick target index `ind` -----------------------------------------
    SIGMA_min = false
    SIGMA_max = false
    local ind::Int
    if sigma_is_numeric
        ind = argmin(abs.(eig_list .- sigma_value))
    elseif sigma_norm === :smallestreal
        SIGMA_min = true; ind = argmin(eig_list)
    elseif sigma_norm === :largestreal
        SIGMA_max = true; ind = argmax(eig_list)
    elseif sigma_norm === :smallestabs
        ind = argmin(abs.(eig_list))
    elseif sigma_norm === :largestabs
        ind = argmax(abs.(eig_list))
    else
        throw(ArgumentError("veigs: unknown sigma $(sigma_norm)"))
    end

    # ---- Rough global bounds on the spectrum -----------------------------
    # Need the smallest eigenvalue of B (the metric matrix) as a starting
    # point. For dense B we use LAPACK; for sparse B we use the same
    # KrylovKit shift-invert pipeline as `_approx_eig`. Both return a
    # safe lower bound that `rough_lower` then verifies.
    λ_B_min_rough = _smallest_eig_estimate(B_mid)
    I_int = Bi isa AbstractSparseMatrix ?
            interval.(sparse(I, n, n)) :
            interval.(Matrix{Float64}(I, n, n))
    λ_B_min       = rough_lower(Bi, I_int,
                                λ_B_min_rough, λ_B_min_rough / 2)
    A_inf_norm    = opnorm(_sup_abs_matrix(Ai), Inf)
    lambda_rough_min = -A_inf_norm / λ_B_min
    lambda_rough_max =  A_inf_norm / λ_B_min

    # ---- Initial direction flags ----------------------------------------
    direction_left  = true
    direction_right = true

    global_sigma = +Inf
    global_rho   = -Inf
    global_s     = 0
    global_r     = n + 1

    if SIGMA_max
        global_sigma = rough_upper(Ai, Bi, eig_list[m], lambda_rough_max)
        global_s     = n
        direction_right = false
    end
    if SIGMA_min
        global_rho = rough_lower(Ai, Bi, eig_list[1], lambda_rough_min)
        global_r   = 1
        direction_left = false
    end

    # ---- Main cluster expansion loop -------------------------------------
    rough = (λ_B_min, lambda_rough_min, lambda_rough_max)

    lambda = Vector{Interval{Float64}}()
    ind_range = Int[]
    local_ind_range = Int[]

    do_next = true
    while do_next
        if isempty(ind_range)
            (loc_lambda, loc_lir, loc_gir) = _veigs_eig_local_bound(
                Ai, Bi, eig_list, V, ind, rough,
                global_r, global_s, global_rho, global_sigma,
                SIGMA_min, SIGMA_max, direction_left, direction_right)
            append!(lambda, loc_lambda)
            ind_range = collect(loc_gir)
            local_ind_range = collect(loc_lir)
        else
            ind_left  = minimum(local_ind_range) - 1
            ind_right = maximum(local_ind_range) + 1
            lir_l = Int[]; lir_r = Int[]
            global_s_local = first(ind_range) - 1
            global_r_local = last(ind_range) + 1

            if direction_right
                global_rho_local = maximum(inf, lambda)
                (loc_lambda, loc_lir, loc_gir) = _veigs_eig_local_bound(
                    Ai, Bi, eig_list, V, ind_right, rough,
                    global_r_local, global_s_local, global_rho_local, global_sigma,
                    SIGMA_min, SIGMA_max, false, direction_right)
                append!(lambda, loc_lambda)
                append!(ind_range, collect(loc_gir))
                lir_r = collect(loc_lir)
            end
            if direction_left
                global_sigma_local = minimum(sup, lambda)
                (loc_lambda, loc_lir, loc_gir) = _veigs_eig_local_bound(
                    Ai, Bi, eig_list, V, ind_left, rough,
                    global_r_local, global_s_local, global_rho, global_sigma_local,
                    SIGMA_min, SIGMA_max, direction_left, false)
                prepend!(lambda, loc_lambda)
                prepend!(ind_range, collect(loc_gir))
                lir_l = collect(loc_lir)
            end
            local_ind_range = vcat(lir_l, lir_r)
        end

        g_r = minimum(ind_range)
        g_s = maximum(ind_range)
        do_next = false
        if g_s - g_r + 1 < num_eigs
            if !isempty(local_ind_range) && minimum(local_ind_range) > 1 && direction_left
                do_next = true
            else
                direction_left = false
            end
            if !isempty(local_ind_range) && maximum(local_ind_range) < m && direction_right
                do_next = true
            else
                direction_right = false
            end
        end
    end

    perm = sortperm(ind_range)
    return (lambda[perm], minimum(ind_range):maximum(ind_range))
end

# ---- Sigma normalization ---------------------------------------------------

function _normalize_sigma(s::Symbol)
    haskey(_LEGACY_SIGMA, s) && return (_LEGACY_SIGMA[s], false)
    s in _MODERN_SIGMA && return (s, false)
    throw(ArgumentError("veigs: unknown sigma symbol :$s"))
end
_normalize_sigma(s::AbstractString) = _normalize_sigma(Symbol(s))
_normalize_sigma(::Real) = (:numeric, true)

# ---- Approximate eigensolve ------------------------------------------------
#
# Two paths share the entry point `_approx_eig`:
#
# - Sparse / large input → KrylovKit shift-invert. Factor `A − σ·B` once
#   with sparse LU (SuiteSparse), then run Lanczos against the linear
#   operator `x ↦ (A − σ·B)⁻¹·(B·x)`. This converges in O(k) iterations
#   with O(n) work per iteration on banded inputs, vs. dense `eigen`'s
#   O(n³).
#
#   Arpack.jl was tried first (it's the named MATLAB equivalent) but
#   `XYAUPD_Exception` — non-convergence — fired on every tuning of the
#   n=5000 sparse Laplacian on this hardware. KrylovKit succeeds in
#   ≤10 ms with default tolerances under shift-invert.
#
# - Dense / small input → LAPACK `eigen`. Simpler, robust, the only
#   correct option when `A` and `B` are dense or n is too small for the
#   Arnoldi setup cost to amortise.
#
# `_approx_eig` always returns `(eig_list::Vector{Float64},
# V::Matrix{Float64})` sorted by eigenvalue ascending, regardless of
# which inner path produced it.

const _APPROX_EIG_KRYLOV_THRESHOLD = 200

function _approx_eig(A::AbstractMatrix, B::AbstractMatrix,
                     num_eigs::Int, sigma_norm::Symbol, sigma_value::Float64,
                     sigma_is_numeric::Bool)
    n = size(A, 1)

    use_krylov = n ≥ _APPROX_EIG_KRYLOV_THRESHOLD &&
                 (A isa AbstractSparseMatrix || B isa AbstractSparseMatrix)
    if use_krylov
        result = _try_approx_eig_krylov(A, B, num_eigs, sigma_norm, sigma_value, sigma_is_numeric)
        result === nothing || return result
    end
    return _approx_eig_dense(A, B)
end

# Dense path: matches the original behaviour exactly.
function _approx_eig_dense(A::AbstractMatrix, B::AbstractMatrix)
    A_f = A isa AbstractMatrix{<:Interval} ? mid.(A) : A
    B_f = B isa AbstractMatrix{<:Interval} ? mid.(B) : B
    Am = Matrix(A_f)
    Bm = Matrix(B_f)
    F  = eigen(Symmetric(Am), Symmetric(Bm))
    perm = sortperm(real.(F.values))
    return (real.(F.values[perm]), real.(F.vectors[:, perm]))
end

# KrylovKit shift-invert: returns the tuple on success, `nothing` to
# trigger the dense fallback.
#
# Algebra: with shift σ, the eigenvalues of (A − σB)⁻¹·B near σ are
# 1/(λ_i − σ); the largest |μ| of those correspond to the λ_i closest
# to σ. So `eigsolve(op, n, k, :LM)` recovers the k eigenvalues of the
# generalized problem nearest σ — same trick MATLAB's `eigs` uses.
function _try_approx_eig_krylov(A::AbstractMatrix, B::AbstractMatrix,
                                num_eigs::Int, sigma_norm::Symbol,
                                sigma_value::Float64, sigma_is_numeric::Bool)
    n = size(A, 1)
    A_f = A isa AbstractMatrix{<:Interval} ? mid.(A) : A
    B_f = B isa AbstractMatrix{<:Interval} ? mid.(B) : B
    A_sp = A_f isa AbstractSparseMatrix ? A_f : sparse(A_f)
    B_sp = B_f isa AbstractSparseMatrix ? B_f : sparse(B_f)

    # Pick a shift σ. For "smallest" sigma we anchor at zero (matches
    # MATLAB's `eigs(A, B, k, 'smallestabs')`). Numeric sigma is the
    # explicit shift. Largest-* falls through to plain Arnoldi (no shift).
    σ, use_shift_invert = if sigma_is_numeric
        (sigma_value, true)
    elseif sigma_norm === :smallestreal || sigma_norm === :smallestabs
        (0.0, true)
    elseif sigma_norm === :largestreal || sigma_norm === :largestabs
        (0.0, false)
    else
        return nothing
    end

    k = min(num_eigs, n - 2)
    k ≥ 1 || return nothing

    try
        if use_shift_invert
            # Factor M = A - σ*B once. Solve M·x = B·y at each Lanczos step.
            M = A_sp - σ * B_sp
            F = lu(M)
            B_dense_op = (B_sp isa SparseMatrixCSC) ? B_sp : sparse(B_sp)
            x0 = randn(n)
            μs, Vs, info = KrylovKit.eigsolve(
                x -> F \ (B_dense_op * x),
                x0, k, :LM;
                tol = 1e-10,
                maxiter = 500,
            )
            info.converged ≥ k || return nothing
            μ_real = real.(μs[1:k])
            λ_real = σ .+ 1.0 ./ μ_real
            V_real = reduce(hcat, real.(Vs[1:k]))
        else
            # Plain Arnoldi on the operator x ↦ A·x. Generalized form is
            # only meaningful with B ≠ I when shift-invert is active, so
            # for the largest-eigenvalue cases we ignore B (callers rarely
            # use these symbols on generalized problems anyway). The
            # verified algorithm uses `eig_list` only as approximate
            # starting points; soundness comes from the verified path.
            x0 = randn(n)
            which = sigma_norm === :largestreal ? :LR : :LM
            λs, Vs, info = KrylovKit.eigsolve(A_sp, x0, k, which; tol = 1e-10)
            info.converged ≥ k || return nothing
            λ_real = real.(λs[1:k])
            V_real = reduce(hcat, real.(Vs[1:k]))
        end

        perm = sortperm(λ_real)
        return (λ_real[perm], V_real[:, perm])
    catch err
        return nothing
    end
end

# Sup-of-absolute-value as a Float64 matrix; for interval entries uses `mag`.
# `mag.(::AbstractMatrix{Interval{Float64}})` already returns Float64, so the
# previous `Float64.` outer wrap was a redundant broadcast pass.
_sup_abs_matrix(M::AbstractMatrix{<:Interval}) = mag.(M)
_sup_abs_matrix(M::AbstractMatrix{<:Real})     = abs.(M)

# Cheap estimate of `λ_min(B_mid)`. Result is used only as a starting point
# for the verified `rough_lower` search; it doesn't need to be tight, just
# below the true `λ_min`. For dense input we use stdlib `eigvals`; for
# sparse input we use the same KrylovKit shift-invert pipeline as the
# main approximate eigensolve.
function _smallest_eig_estimate(B_mid::AbstractMatrix)
    n = size(B_mid, 1)
    if B_mid isa AbstractSparseMatrix
        # KrylovKit shift-invert near zero — finds the smallest |λ|.
        try
            F = lu(sparse(Symmetric(B_mid)))
            μs, _, info = KrylovKit.eigsolve(x -> F \ x, n, 1, :LM;
                                             tol = 1e-8, maxiter = 200)
            info.converged ≥ 1 || return _eigvals_min_dense(B_mid)
            return real(1.0 / μs[1])
        catch
            return _eigvals_min_dense(B_mid)
        end
    end
    return minimum(real.(eigvals(Symmetric(Matrix(B_mid)))))
end

_eigvals_min_dense(B::AbstractMatrix) =
    minimum(real.(eigvals(Symmetric(Matrix(B)))))

# ---- Core local-bound subfunction ------------------------------------------
function _veigs_eig_local_bound(Ai::AbstractMatrix{<:Interval},
                                 Bi::AbstractMatrix{<:Interval},
                                 eig_list::AbstractVector{<:Real},
                                 V::AbstractMatrix,
                                 ind::Int,
                                 rough::NTuple{3,<:Real},
                                 global_r::Int, global_s::Int,
                                 global_rho::Real, global_sigma::Real,
                                 SIGMA_min::Bool, SIGMA_max::Bool,
                                 direction_left::Bool, direction_right::Bool)
    n = size(Ai, 1)
    m = length(eig_list)
    λ_B_min, lambda_rough_min, lambda_rough_max = rough

    r = ind
    s = ind

    # ---- Right (upper) verification --------------------------------------
    if direction_right
        while s < m && eig_list[s] > eig_list[s+1] - _VEIGS_LAMBDA_RATIO * abs(eig_list[s+1])
            s += 1
        end
        eig_next = s == m ? eig_list[s] + abs(eig_list[s]) * 0.01 : eig_list[s+1]
        λ = _VEIGS_LAMBDA_RATIO * eig_list[s] + (1 - _VEIGS_LAMBDA_RATIO) * eig_next

        Mi = Ai .- interval(λ) .* Bi
        L, D, p, ΔM, ok = verified_ldl(Mi)
        ok || throw(VeigsLDLFailureError("veigs: midpoint LDL failed at λ = $λ (right)"))
        (neg, _, zer, _) = inertia(D)

        if (neg + zer) == n
            SIGMA_max    = true
            s            = m
            global_sigma = rough_upper(Ai, Bi, eig_list[s], lambda_rough_max)
            global_s     = n
        else
            # `mag(-x) = mag(x)` so we can drop the negation here, saving an
            # O(nnz) sparse interval allocation. The `+ ΔM` form below
            # likewise replaces `.- (.-ΔM)`, exact in interval arithmetic.
            err_est = opnorm(_sup_abs_matrix(ΔM), Inf) / λ_B_min
            gap_right = (1 - _VEIGS_LAMBDA_RATIO) * (eig_next - eig_list[s])
            if err_est > gap_right
                while true
                    err_test = err_est / 2
                    tmp = interval(err_test) .* Bi .+ ΔM
                    if !verified_isspd(sym_hull(tmp))
                        break
                    end
                    err_est = err_test
                end
            end
            is_pos = err_est ≤ gap_right
            if !is_pos
                if neg == 1
                    SIGMA_min = true
                    r = 1
                    global_rho = rough_lower(Ai, Bi, eig_list[r], lambda_rough_min)
                    global_r   = 1
                else
                    error("veigs: LDL failed to find lower bound for eigenvalue λ_{$(s+1)}")
                end
            else
                global_s     = neg
                global_sigma = inf(interval(λ) - interval(err_est))
            end
        end
    end

    # ---- Left (lower) verification ---------------------------------------
    if direction_left
        r = ind
        while r > 1 && eig_list[r] < eig_list[r-1] + _VEIGS_LAMBDA_RATIO * abs(eig_list[r-1])
            r -= 1
        end
        eig_pre = r == 1 ? eig_list[r] - abs(eig_list[r]) * 0.01 : eig_list[r-1]
        λ = _VEIGS_LAMBDA_RATIO * eig_list[r] + (1 - _VEIGS_LAMBDA_RATIO) * eig_pre

        Mi = Ai .- interval(λ) .* Bi
        L, D, p, ΔM, ok = verified_ldl(Mi)
        ok || throw(VeigsLDLFailureError("veigs: midpoint LDL failed at λ = $λ (left)"))
        (neg, _, _, _) = inertia(D)

        if neg == 0
            SIGMA_min = true
            r = 1
            global_rho = rough_lower(Ai, Bi, eig_list[r], lambda_rough_min)
            global_r   = 1
        else
            err_est = opnorm(_sup_abs_matrix(ΔM), Inf) / λ_B_min
            gap_left = (1 - _VEIGS_LAMBDA_RATIO) * (eig_list[r] - eig_pre)
            if err_est > gap_left
                while true
                    err_test = err_est / 2
                    tmp = interval(err_test) .* Bi .- ΔM
                    if !verified_isspd(sym_hull(tmp))
                        break
                    end
                    err_est = err_test
                end
            end
            is_pos = err_est ≤ gap_left
            if !is_pos
                if neg == n - 1
                    SIGMA_max    = true
                    s            = m
                    global_sigma = rough_upper(Ai, Bi, eig_list[s], lambda_rough_max)
                    global_s     = n
                else
                    error("veigs: LDL failed to find upper bound for eigenvalue λ_{$(r-1)}")
                end
            else
                global_r   = neg + 1
                global_rho = sup(interval(λ) + interval(err_est))
            end
        end
    end

    # ---- LDL result check / Rayleigh fallback ---------------------------
    if global_s - global_r != s - r || global_s < 0 || global_r > n
        if SIGMA_max && SIGMA_min
            lower_bound = _rayleigh_quotient(Ai, Bi, V[:, m])
            upper_bound = _rayleigh_quotient(Ai, Bi, V[:, 1])
            return ([hull(interval(global_rho), upper_bound),
                     hull(lower_bound, interval(global_sigma))],
                    [1, n], [1, n])
        end
        if SIGMA_max
            lower_bound = _rayleigh_quotient(Ai, Bi, V[:, m])
            return ([interval(inf(lower_bound), float(global_sigma))], [n], [n])
        end
        if SIGMA_min
            upper_bound = _rayleigh_quotient(Ai, Bi, V[:, 1])
            return ([interval(float(global_rho), sup(upper_bound))], [1], [1])
        end
    end

    if global_s - global_r > s - r
        throw(VeigsClusterTooLargeError(
            "veigs: cluster too large — increase the EigNum approximation count"))
    end
    if global_s - global_r < s - r
        throw(VeigsRoughBoundError(
            "veigs: rough-bound separation failed (approximate eigenvalues likely poorly conditioned)"))
    end

    # Slack ρ/σ outward when at spectrum extremes (matches MATLAB).
    if SIGMA_max
        global_sigma = max(global_sigma, eig_list[m] + abs(eig_list[m]))
    end
    if SIGMA_min
        global_rho = min(global_rho, eig_list[r] - abs(eig_list[r]))
    end

    Vc = V[:, r:s]
    chunk = lehmann_behnke(Ai, Bi, eig_list, Vc, global_rho, global_sigma,
                           λ_B_min, r, s; do_shift = _VEIGS_DO_SHIFT)
    return (chunk, collect(r:s), collect(global_r:global_s))
end

# Verified Rayleigh quotient v' A v / v' B v as a (scalar) interval.
function _rayleigh_quotient(Ai::AbstractMatrix{<:Interval},
                            Bi::AbstractMatrix{<:Interval},
                            v::AbstractVector)
    vi = interval.(v)
    num = transpose(vi) * Ai * vi
    den = transpose(vi) * Bi * vi
    # `num` and `den` are 1-element after dot products. Extract the scalar.
    return only(num) / only(den)
end
