# src/veig.jl
#
# Port of MATLAB `veig.m` (Liu/Yanagisawa). Yamamoto-style verified solver
# for the small-to-medium dense generalized symmetric eigenproblem
# `A x = λ B x`, with `B` (interval) positive definite.
#
# CONTRACT (soundness):
#   `veig(A, B, ind)` returns `(eig_bounds, ind_range)` where every concrete
#   matrix pair (A_c, B_c) inside the interval enclosure of (A, B) has its
#   eigenvalue λ_i (for each `i ∈ ind_range`) inside `eig_bounds[i - first(ind_range) + 1]`.
#   The bounds may *cluster*: if two eigenvalues are too close to certify
#   separately, both indices are mapped to the same wider interval.
#
# Algorithm (from veig.m):
#   1. Approximate eigenvalues `eig_list = eigvals(mid(A), mid(B))` (sorted).
#   2. For each requested index, search rough lower/upper test points by
#      counting LDL inertia of `mid(A) - λ_test mid(B)`. Doubling-decimal
#      schedule `λ ± max(eps(λ_test), |λ_test|) * 10^(k-10)` capped at k=10.
#   3. Verified lower bound: factor `mid(A - λ_test_low B)` by Bunch-Kaufman,
#      compute interval residual `ΔM = (A - λ_test_low B) - PLDLᵀPᵀ`, then
#      double `err_est` until both `err_est B + ΔM` and `err_est B - ΔM` are
#      verified PD. Cap k=53 (matches MATLAB).
#   4. Symmetric upper-bound branch.
#   5. `eig_bound = hull(λ_lower, λ_upper)`; tag every index in
#      `(neg_num_low+1) : neg_zero_num_upper` with this bound.
#
# Per Q-004, the n > 100 cap is kept for parity with MATLAB.

using LinearAlgebra
using SparseArrays
using IntervalArithmetic: Interval, interval, mid, hull, mag

const _VEIG_MAX_ROUGH_K = 30        # rough-bound iter cap. MATLAB has 10
                                    # (veig.m lines 114, 145); we raise to 30
                                    # because veig is called on small Lehmann
                                    # subproblems whose eigenvalue magnitudes
                                    # can be far off 1.0, and the relative
                                    # step schedule abs(eig_test)*10^(k-10)
                                    # needs more headroom in those cases.
                                    # See D-017.
const _VEIG_MAX_VERIFIED_K = 53     # verified-bound iter cap (MATLAB line 182, 218)
const _VEIG_SIZE_CAP = 100          # MATLAB hard cap (Q-004; keep for parity)

"""
    veig(A, B[, ind]) -> (eig_bounds, ind_range)

Verified bounds for the eigenvalues of the symmetric generalized eigenproblem
`A x = λ B x`, indexed by ascending eigenvalue order. `A` and `B` may be real
or interval matrices. `B` must be positive definite (interval PD allowed).

`ind` selects the eigenvalue indices to bound (default: all of them). Pass a
single integer or any iterable yielding an integer range.

Returns:
- `eig_bounds::Vector{Interval{Float64}}` — verified enclosures.
- `ind_range::UnitRange{Int}` — the eigenvalue indices actually covered.
  When two eigenvalues cannot be separated, both indices map to the same
  (wider) bound; `ind_range` may extend slightly beyond the requested `ind`.

Throws:
- `VeigsSizeError` if `A`, `B` are not the same square shape, are not
  symmetric in midpoint, or `n > 100` (Q-004).
- `VeigsLDLFailureError` if a midpoint LDL factorization fails inside the
  bound-search loops.
- `ErrorException` if the rough-bound or verified-bound iteration caps are
  hit (matches MATLAB error sites).
"""
function veig(A::AbstractMatrix, B::AbstractMatrix, ind = 1:size(A, 1))
    n = size(A, 1)
    n == size(A, 2) ||
        throw(VeigsSizeError("veig: A must be square"))
    size(B) == size(A) ||
        throw(VeigsSizeError("veig: A and B must have the same shape"))
    n ≤ _VEIG_SIZE_CAP ||
        throw(VeigsSizeError("veig: matrix too large (n=$n > $_VEIG_SIZE_CAP); use the cluster-aware veigs driver"))

    # Promote to interval (idempotent if already interval)
    Ai = A isa AbstractMatrix{<:Interval} ? A : interval.(Matrix(A))
    Bi = B isa AbstractMatrix{<:Interval} ? B : interval.(Matrix(B))
    A_mid = mid.(Ai)
    B_mid = mid.(Bi)

    # Midpoint symmetry — MATLAB checks `isequal(A, A')` strictly. Allow tiny
    # floating-point asymmetry but reject anything beyond eps-relative.
    _check_sym_midpoint(A_mid, "A")
    _check_sym_midpoint(B_mid, "B")

    # Index range
    isempty(ind) && throw(ArgumentError("veig: ind must be non-empty"))
    min_ind, max_ind = extrema(ind)
    (min_ind ≥ 1 && max_ind ≤ n) ||
        throw(ArgumentError("veig: ind must satisfy 1 ≤ ind ≤ $n, got [$min_ind, $max_ind]"))

    # Approximate eigenvalues, sorted ascending. Symmetric problem.
    eig_list = sort!(real.(eigvals(Symmetric(A_mid), Symmetric(B_mid))))
    min_eig_B = minimum(real.(eigvals(Symmetric(B_mid))))

    eig_bounds = Vector{Interval{Float64}}()
    local_min_ind = min_ind     # may grow downward as clusters absorb earlier indices

    index = min_ind
    while index ≤ max_ind
        eig_test = eig_list[index]

        # ---- Rough lower bound ---------------------------------------------
        if index == 1
            neg_num_low     = 0
            lambda_test_low = eig_test
        else
            (neg_num_low, lambda_test_low) =
                _rough_lower_search(A_mid, B_mid, eig_test, index)
        end

        # ---- Rough upper bound ---------------------------------------------
        if index == n
            neg_zero_num_upper = n
            lambda_test_upper  = eig_test
        else
            (neg_zero_num_upper, lambda_test_upper) =
                _rough_upper_search(A_mid, B_mid, eig_test, index)
        end

        # ---- Verified lower bound ------------------------------------------
        lambda_lower = _verified_one_sided_bound(
            Ai, Bi, lambda_test_low, eig_test, min_eig_B; side = :lower)

        # ---- Verified upper bound ------------------------------------------
        lambda_upper = _verified_one_sided_bound(
            Ai, Bi, lambda_test_upper, eig_test, min_eig_B; side = :upper)

        # ---- Store bound ---------------------------------------------------
        eig_bound = hull(interval(lambda_lower), interval(lambda_upper))
        index_range = (neg_num_low + 1):neg_zero_num_upper

        # If the cluster reaches below our current local_min_ind, extend
        # downward (matches MATLAB lines 241-243).
        if neg_num_low + 1 < local_min_ind
            local_min_ind = neg_num_low + 1
        end

        # Resize and broadcast eig_bound across `index_range`. Indices in the
        # returned vector are 1-based starting at `local_min_ind`.
        for k in index_range
            slot = k - local_min_ind + 1
            while length(eig_bounds) < slot
                # Placeholder; always overwritten before return because every
                # slot in the final `ind_range` is touched by some cluster.
                push!(eig_bounds, interval(0.0, 0.0))
            end
            eig_bounds[slot] = eig_bound
        end

        index = neg_zero_num_upper + 1
    end

    ind_range = local_min_ind:(index - 1)
    return (eig_bounds, ind_range)
end

# ---- Helpers ---------------------------------------------------------------

function _check_sym_midpoint(M::AbstractMatrix, name::String)
    tol = eps(eltype(M)) * (1 + maximum(abs, M))
    err = maximum(abs, M - transpose(M))
    err ≤ tol ||
        throw(VeigsSizeError("veig: matrix $name is not symmetric (max |M - M'| = $err > tol $tol)"))
    return nothing
end

# Float-only midpoint LDL + inertia. Used inside the rough-bound search.
# Returns (neg, pos, zer, F::Bool) and a `failed::Bool` flag.
# `failed=true` only for hard LAPACK failures (info < 0). Singular pivots
# (info > 0) still yield valid factors with a zero block in D, which
# `inertia` correctly counts as a zero eigenvalue.
function _midpoint_inertia(M_mid::AbstractMatrix)
    Mm = (M_mid + transpose(M_mid)) / 2     # defensive symmetrization
    F = bunchkaufman(Symmetric(Mm, :L); check = false)
    F.info < 0 && return (0, 0, 0, true, true)   # true LAPACK failure
    D = Matrix(F.D)
    (neg, pos, zer, F_unc) = inertia(D)
    return (neg, pos, zer, F_unc, false)
end

function _rough_lower_search(A_mid::AbstractMatrix, B_mid::AbstractMatrix,
                             eig_test::Real, index::Int)
    lambda_test = float(eig_test)
    for k in 0:_VEIG_MAX_ROUGH_K
        E = A_mid .- lambda_test .* B_mid
        (neg, _, _, F_unc, failed) = _midpoint_inertia(E)
        if failed || F_unc || neg ≥ index
            # Step further down.
            step = max(eps(eig_test), abs(eig_test)) * 10.0^(k - 10)
            lambda_test = float(eig_test) - step
        else
            return (neg, lambda_test)
        end
    end
    error("veig: rough lower bound search exceeded iteration limit at index $index")
end

function _rough_upper_search(A_mid::AbstractMatrix, B_mid::AbstractMatrix,
                             eig_test::Real, index::Int)
    lambda_test = float(eig_test)
    for k in 0:_VEIG_MAX_ROUGH_K
        E = A_mid .- lambda_test .* B_mid
        (neg, _, zer, F_unc, failed) = _midpoint_inertia(E)
        if failed || F_unc || (neg + zer) < index
            step = max(eps(eig_test), abs(eig_test)) * 10.0^(k - 10)
            lambda_test = float(eig_test) + step
        else
            return (neg + zer, lambda_test)
        end
    end
    error("veig: rough upper bound search exceeded iteration limit at index $index")
end

# Verified-bound branch (lower or upper). Returns the certified scalar bound.
function _verified_one_sided_bound(Ai::AbstractMatrix{<:Interval},
                                   Bi::AbstractMatrix{<:Interval},
                                   lambda_test::Real, eig_test::Real,
                                   min_eig_B::Real; side::Symbol)
    side === :lower || side === :upper ||
        throw(ArgumentError("side must be :lower or :upper"))

    # Construct the interval matrix M = A - lambda_test * B and factor it.
    λi = interval(lambda_test)
    Mi = Ai .- λi .* Bi
    L, D, p, ΔM, ok = verified_ldl(Mi)
    ok || throw(VeigsLDLFailureError(
        "veig: midpoint LDL of (A - λ B) failed at λ = $lambda_test ($side branch)"))

    # diff_m_inf = ‖ΔM‖_∞ via mag (interval-aware sup of |·|).
    diff_m_inf = opnorm(mag.(ΔM), Inf)

    if diff_m_inf == 0
        return float(lambda_test)
    end

    err_est = max(eps(eig_test), diff_m_inf / min_eig_B)
    is_positive = false
    for _ in 1:_VEIG_MAX_VERIFIED_K
        err_est *= 2
        err_est_i = interval(err_est)
        # Both directions: err_est*B ± ΔM must be PD.
        tmp_plus  = err_est_i .* Bi .+ ΔM
        tmp_minus = err_est_i .* Bi .- ΔM
        if verified_isspd(sym_hull(tmp_plus)) && verified_isspd(sym_hull(tmp_minus))
            is_positive = true
            break
        end
    end
    is_positive || error(
        "veig: failed to certify err_est in verified $(side) bound at λ = $lambda_test")

    return side === :lower ? float(lambda_test) - err_est : float(lambda_test) + err_est
end
