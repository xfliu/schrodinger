module Veigs

# ============================================================================
# Veigs — verified eigenvalue bounds for the symmetric generalized eigenproblem
# A x = λ B x.
#
# Public API mirrors the MATLAB `veigs` toolbox:
#   veigs(A, B[, k][, sigma]) -> (lambda, ind_range)
#   veig(A, B[, ind])         -> (eig_bounds, ind_range)
#
# Internal building blocks (also exported for cross-package reuse and testing):
#   inertia, verified_isspd, sym_hull
#
# See DECISIONS.md (next to this file) for the active design choices.
# ============================================================================

using LinearAlgebra
using SparseArrays
using IntervalArithmetic: Interval, interval, mid, hull, radius, mag,
                           inf, sup

# ---- Error hierarchy (D-008) -----------------------------------------------
abstract type VeigsError <: Exception end

struct VeigsLDLFailureError      <: VeigsError; msg::String; end
struct VeigsClusterTooLargeError <: VeigsError; msg::String; end
struct VeigsLehmannBehnkeError   <: VeigsError; msg::String; end
struct VeigsRoughBoundError      <: VeigsError; msg::String; end
struct VeigsSizeError            <: VeigsError; msg::String; end

Base.showerror(io::IO, e::VeigsError) = print(io, typeof(e), ": ", e.msg)

# ---- Tiny helpers ----------------------------------------------------------
"""
    sym_hull(M; assume_symmetric=false) -> AbstractMatrix
    sym_hull!(H::AbstractMatrix, M::AbstractMatrix) -> H

Symmetric hull: each entry of the result is `hull(M[i,j], M[j,i])`.
Mirrors MATLAB `hull(M, M')`. Used wherever the MATLAB code defends
against asymmetric round-off in interval matrices before calling
`verified_isspd`.

Implemented as a tight upper-triangle loop (only `n(n+1)/2` calls to
`hull`, vs `n²` for the broadcast form `hull.(M, transpose(M))`). At
n=5000 the broadcast form was ~21 % of total `veigs` wall time on the
profile of 2026-05-01 — most of that overhead is the broadcast
machinery and `transpose` indexing, not the per-element `hull`.

`assume_symmetric=true` skips the work entirely and returns `M` (or its
guaranteed-symmetric equivalent). Pass it only when the caller knows
the matrix is already symmetric in interval arithmetic — typically when
`M` was constructed by an expression like `A − λ·B` of two symmetric
matrices, where the broadcast preserves entry-wise equality
`M[i,j] == M[j,i]`. Use the default (`false`) for matrices built from
non-symmetric expressions like `V'·A·nV − nV'·(B·nV − A·V) + Err`,
where interval round-off can introduce asymmetry.
"""
function sym_hull(M::AbstractMatrix; assume_symmetric::Bool = false)
    assume_symmetric && return M
    n, m = size(M)
    n == m || throw(DimensionMismatch("sym_hull: matrix must be square, got $(n)×$(m)"))
    H = similar(M)
    return sym_hull!(H, M)
end

function sym_hull!(H::AbstractMatrix, M::AbstractMatrix)
    n = size(M, 1)
    @inbounds for j in 1:n
        H[j, j] = hull(M[j, j], M[j, j])
        for i in 1:(j - 1)
            h = hull(M[i, j], M[j, i])
            H[i, j] = h
            H[j, i] = h
        end
    end
    return H
end

"""
    _check_finite_midpoint(M)

Throws `DomainError` if any entry of `mid.(M)` is `Inf` or `NaN`.
Used by `verified_isspd` (D-015) and other entry points.
"""
function _check_finite_midpoint(M::AbstractMatrix)
    Mm = mid.(M)
    if any(!isfinite, Mm)
        throw(DomainError(M, "matrix contains Inf or NaN entries"))
    end
    return nothing
end

"""
    _midrad_floats(M) -> (mid::AbstractMatrix{Float64}, rad)

Split `M` into a float midpoint matrix and either a float radius matrix or
`nothing` (meaning "radius identically zero", for exact float input).
Sparsity-preserving on `SparseMatrixCSC` input. Used by primitives that
follow DESIGN.md §1: operate on `(mid, rad)` pairs in float without ever
materialising an n×n interval matrix.
"""
function _midrad_floats(M::AbstractMatrix)
    if M isa AbstractMatrix{<:Interval}
        return Float64.(mid.(M)), Float64.(radius.(M))
    elseif M isa AbstractSparseMatrix
        # Avoid `Float64.(::SparseMatrixCSC)` triggering an eltype copy when
        # already Float64; otherwise convert.
        return (eltype(M) === Float64 ? M : sparse(Float64.(M))), nothing
    else
        return (eltype(M) === Float64 ? M : Matrix{Float64}(M)), nothing
    end
end

# ---- Submodule includes (each routine is its own file) ---------------------
include("inertia.jl")
include("verified_ldl.jl")
include("verified_isspd.jl")    # depends on verified_ldl + inertia
include("rough_bounds.jl")
include("veig.jl")
include("lehmann_behnke.jl")
include("veigs_driver.jl")     # iter 7 — `veigs` function lives here

export veig, veigs,
       inertia, verified_isspd, verified_ldl, rough_lower, rough_upper,
       lehmann_behnke, sym_hull, DEFAULT_ISSPD_METHOD,
       VeigsError, VeigsLDLFailureError, VeigsClusterTooLargeError,
       VeigsLehmannBehnkeError, VeigsRoughBoundError, VeigsSizeError

end # module
