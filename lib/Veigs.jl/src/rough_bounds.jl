# src/rough_bounds.jl
#
# Port of MATLAB `veigs_RoughLower` / `veigs_RoughUpper` (subfunctions in
# `veigs.m`, lines ~518-570). Search exponentially outward from a seed
# eigenvalue toward a floor/ceiling, calling `verified_isspd` on each
# candidate `A - λB` (or `λB - A`) until we find a verified-PD shift.
#
# CONTRACT (soundness):
#   - On success, `rough_lower(A, B, λ_seed, λ_floor) == λ` certifies that
#     `A - λB` is positive definite, i.e. every eigenvalue of `A x = μ B x`
#     satisfies `μ > λ`. Equivalently, `λ` is a verified lower bound on the
#     entire spectrum.
#   - On success, `rough_upper(A, B, λ_seed, λ_ceiling) == λ` certifies that
#     `λB - A` is positive definite, so every spectral eigenvalue is `< λ`.
#   - On failure (no λ in [floor, seed] or [seed, ceiling] verified PD),
#     throws `VeigsRoughBoundError`.
#
# Schedule (matches MATLAB exactly per guide.md §4 iter 4):
#   delta = max(eps, |seed - bound_extreme| / 2^32)
#   k = 1, 2, 3, ...
#   λ = seed ± delta * 2^(k-1)
# The first iteration (k=1) always runs even if `seed` already crosses the
# bound — matches the MATLAB `(lambda >= eig_test_min - eps) || (k == 1)`
# loop guard.

using LinearAlgebra
using IntervalArithmetic: Interval, interval

# Hard cap on the number of doublings. MATLAB has no explicit cap, but the
# schedule reaches λ = bound at k=33 by construction (delta * 2^32 ==
# |seed - bound|), so any value beyond ~40 is purely defensive.
const _ROUGH_BOUNDS_MAX_K = 64

"""
    rough_lower(A, B, λ_seed, λ_floor) -> Float64

Search downward from `λ_seed` toward `λ_floor` (with `λ_floor < λ_seed`)
for a value `λ` such that `verified_isspd(A - λB)` returns `true`. Doubling
schedule: `λ_k = λ_seed - delta * 2^(k-1)` with `delta = max(eps, (λ_seed - λ_floor)/2^32)`.

Returns the first `λ` at which the test succeeds.

Throws `VeigsRoughBoundError` if no `λ` in `[λ_floor, λ_seed]` certifies PD.
"""
function rough_lower(A::AbstractMatrix, B::AbstractMatrix,
                     λ_seed::Real, λ_floor::Real)
    return _rough_search(A, B, λ_seed, λ_floor; direction = :lower)
end

"""
    rough_upper(A, B, λ_seed, λ_ceiling) -> Float64

Search upward from `λ_seed` toward `λ_ceiling` (with `λ_ceiling > λ_seed`)
for a value `λ` such that `verified_isspd(λB - A)` returns `true`. Schedule
mirrors `rough_lower`.

Throws `VeigsRoughBoundError` if no `λ` in `[λ_seed, λ_ceiling]` certifies PD.
"""
function rough_upper(A::AbstractMatrix, B::AbstractMatrix,
                     λ_seed::Real, λ_ceiling::Real)
    return _rough_search(A, B, λ_seed, λ_ceiling; direction = :upper)
end

# Shared core. `direction = :lower` searches down toward the floor and tests
# `A - λB` for PD; `direction = :upper` searches up toward the ceiling and
# tests `λB - A` for PD.
function _rough_search(A::AbstractMatrix, B::AbstractMatrix,
                       λ_seed::Real, λ_extreme::Real;
                       direction::Symbol)
    direction === :lower || direction === :upper ||
        throw(ArgumentError("direction must be :lower or :upper"))
    size(A) == size(B) ||
        throw(VeigsSizeError("rough_bounds: A and B must have the same shape"))
    size(A, 1) == size(A, 2) ||
        throw(VeigsSizeError("rough_bounds: A must be square"))

    if direction === :lower
        λ_extreme < λ_seed ||
            throw(ArgumentError("rough_lower: λ_floor ($λ_extreme) must be < λ_seed ($λ_seed)"))
    else
        λ_extreme > λ_seed ||
            throw(ArgumentError("rough_upper: λ_ceiling ($λ_extreme) must be > λ_seed ($λ_seed)"))
    end

    # Split A, B into (mid::Float64, rad::Union{Float64-matrix, Nothing}) once.
    # `rad === nothing` means the input was an exact float (zero radius); we
    # skip the input-radius term in that case (DESIGN.md §1: never construct
    # an n×n interval matrix here).
    A_mid, A_rad = _midrad_floats(A)
    B_mid, B_rad = _midrad_floats(B)

    delta = max(eps(), abs(λ_seed - λ_extreme) / 2.0^32)
    u = eps(Float64) / 2
    rfac = u / (1 - u)
    n = size(A_mid, 1)

    # Precompute scalar p-norm bounds on the float A_mid, B_mid (and any
    # input radii). These don't change with λ; reused across all loop iters.
    nrm1_A    = opnorm(A_mid, 1)
    nrmI_A    = opnorm(A_mid, Inf)
    nrm1_B    = opnorm(B_mid, 1)
    nrmI_B    = opnorm(B_mid, Inf)
    nrm1_Arad = A_rad === nothing ? 0.0 : opnorm(A_rad, 1)
    nrmI_Arad = A_rad === nothing ? 0.0 : opnorm(A_rad, Inf)
    nrm1_Brad = B_rad === nothing ? 0.0 : opnorm(B_rad, 1)
    nrmI_Brad = B_rad === nothing ? 0.0 : opnorm(B_rad, Inf)

    for k in 1:_ROUGH_BOUNDS_MAX_K
        step = delta * 2.0^(k - 1)
        λ = direction === :lower ? λ_seed - step : λ_seed + step

        # First iteration always runs (matches MATLAB `|| (k == 1)` clause).
        # Subsequent iterations exit if λ has crossed the bound.
        if k > 1
            crossed = direction === :lower ?
                (λ < λ_extreme - eps()) : (λ > λ_extreme + eps())
            crossed && break
        end

        # Form M_mid in float without ever materialising an interval matrix.
        # Both directions produce a *symmetric* float matrix entry-wise
        # because A_mid, B_mid are symmetric and λ is scalar.
        M_mid = direction === :lower ? A_mid .- λ .* B_mid : λ .* B_mid .- A_mid

        # Scalar p-norm bounds on the radius of M_mid as a representation of
        # the true mathematical residual. Three contributions:
        #   (a) input radius from A     (nrm*_Arad)
        #   (b) input radius from B     (|λ|·nrm*_Brad)
        #   (c) Wilkinson rounding on the float computation, bounded by
        #       u/(1-u)·(‖M_mid‖_p + |λ|·‖B_mid‖_p)
        # Sub-additivity of the norms gives an upper bound for the full
        # radius without ever forming a per-entry M_rad matrix (DESIGN.md §1).
        absλ    = abs(λ)
        nrm1_M  = opnorm(M_mid, 1)
        nrmI_M  = opnorm(M_mid, Inf)
        r1_bnd   = rfac * (nrm1_M + absλ * nrm1_B) + nrm1_Arad + absλ * nrm1_Brad
        rinf_bnd = rfac * (nrmI_M + absλ * nrmI_B) + nrmI_Arad + absλ * nrmI_Brad

        if _verified_isspd_rump2006_midrad(M_mid, (r1_bnd, rinf_bnd), n)
            return λ
        end
    end

    throw(VeigsRoughBoundError(
        string("rough_", direction === :lower ? "lower" : "upper",
               ": no verified-PD λ in ",
               direction === :lower ? "[$λ_extreme, $λ_seed]" : "[$λ_seed, $λ_extreme]")))
end
