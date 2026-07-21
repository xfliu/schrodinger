# src/inertia.jl
#
# Port of MATLAB GetInertia.m (Xuefeng Liu 2011, Yuuka Yanagisawa 2025).
#
# Counts (neg, pos, zero) eigenvalues of a block-diagonal symmetric matrix
# composed of 1×1 and 2×2 blocks (the `D` factor of a Bunch-Kaufman LDLᵀ
# factorization). Uses interval arithmetic on each block's determinant /
# trace so that signs are certified, not merely "computed".
#
# CONTRACT (soundness):
#   - When the returned `F == false`, the (neg, pos, zero) counts are
#     mathematically certified to be the inertia of `D`.
#   - When `F == true`, the routine could not certify the inertia of at
#     least one block (e.g., a 2×2 block with determinant interval
#     straddling zero). The (neg, pos, zero) values are partial and should
#     not be relied on.
#
# See plan §7.1 for the corner-case taxonomy.

using IntervalArithmetic: Interval, interval, inf, sup

"""
    inertia(D; abort_on_err=false) -> (neg::Int, pos::Int, zero::Int, F::Bool)

Compute the inertia of a block-diagonal symmetric matrix `D` made of 1×1
and 2×2 diagonal blocks. `D` may be a real matrix or an interval matrix;
real input is wrapped into intervals internally.

Returns the number of negative, positive, and zero eigenvalues, plus a
failure flag `F`. When `F == false` the counts are certified.

If `abort_on_err == true`, raises `VeigsLDLFailureError` instead of
returning `F = true` on undecidable blocks.

# Examples
```julia
julia> inertia(Diagonal([1.0, -2.0, 0.0]))
(1, 1, 1, false)

julia> inertia([0.0 1.0; 1.0 0.0])      # forced 2×2 block, mixed signs
(1, 1, 0, false)
```
"""
function inertia(D::Diagonal; abort_on_err::Bool = false)
    # Fast path: D is diagonal (no 2×2 blocks possible). Iterate the
    # diagonal directly and avoid the O(n²) `interval.(D)` broadcast that
    # the generic path does for a dense `Matrix{Float64}` D.
    diagD = D.diag
    n = length(diagD)
    neg = 0; pos = 0; zer = 0
    if eltype(diagD) <: Interval
        @inbounds for d in diagD
            if _is_strictly_pos(d)
                pos += 1
            elseif _is_strictly_neg(d)
                neg += 1
            else
                zer += 1
            end
        end
    else
        @inbounds for d in diagD
            if d > 0
                pos += 1
            elseif d < 0
                neg += 1
            else
                zer += 1
            end
        end
    end
    return (neg, pos, zer, false)
end

function inertia(D::AbstractMatrix; abort_on_err::Bool = false)
    n = size(D, 1)
    size(D, 2) == n || throw(VeigsSizeError("inertia: D must be square"))

    # Promote to intervals once at the boundary; everything below is interval.
    Di = D isa AbstractMatrix{<:Interval} ? D : interval.(D)

    neg = 0; pos = 0; zer = 0
    i = 1
    while i ≤ n - 1
        off = Di[i, i+1]
        # ---- 1×1 block ------------------------------------------------------
        if _is_zero(off)
            d = Di[i, i]
            if _is_strictly_pos(d)
                pos += 1
            elseif _is_strictly_neg(d)
                neg += 1
            else
                zer += 1                       # certified zero (or surrounded by 0)
            end
            i += 1
            continue
        end

        # ---- 2×2 block ------------------------------------------------------
        # det = D[i,i]*D[i+1,i+1] - D[i,i+1]*D[i+1,i]
        det_block = Di[i, i] * Di[i+1, i+1] - Di[i, i+1] * Di[i+1, i]

        if sup(det_block) < 0
            # Eigenvalues have different signs -> one positive, one negative.
            pos += 1; neg += 1
            i += 2
            continue
        end

        if inf(det_block) > 0
            # Same sign and non-zero. Trace decides which sign.
            tr_block = Di[i, i] + Di[i+1, i+1]
            if inf(tr_block) > 0
                pos += 2; i += 2; continue
            elseif sup(tr_block) < 0
                neg += 2; i += 2; continue
            end
        end

        # ---- Undecidable ----------------------------------------------------
        if abort_on_err
            throw(VeigsLDLFailureError(
                "inertia: cannot certify sign at block (i=$i): " *
                "det ∈ [$(inf(det_block)), $(sup(det_block))]"))
        else
            return (neg, pos, zer, true)
        end
    end

    # Trailing 1×1 block when n is odd (or when the loop finishes at i==n).
    if i == n
        d = Di[i, i]
        if _is_strictly_pos(d)
            pos += 1
        elseif _is_strictly_neg(d)
            neg += 1
        else
            zer += 1
        end
    end

    # Backfill: anything not counted is zero (matches MATLAB final line).
    zer = n - pos - neg
    return (neg, pos, zer, false)
end

# ---- Tiny sign predicates over intervals ------------------------------------
@inline _is_zero(x::Interval)        = inf(x) == 0 && sup(x) == 0
@inline _is_strictly_pos(x::Interval) = inf(x) > 0
@inline _is_strictly_neg(x::Interval) = sup(x) < 0
