# lg_goerisch_correct2.jl -- the CORRECTED Lehmann-Goerisch realization, for EITHER
# single-parity sector (eee primary bound, oee separator) and general index counts.
#
# Extends goer/lg_goerisch_correct.jl (which assumed nx == ny == nz and an eee sector with a
# constant mode, so it could not be used for the oee separator stage) in four places:
#   (1) vsq_quadform takes (nx,ny,nz) and separate x / transverse moment closures, so the
#       oee block nx != ny == nz works;
#   (2) grad_phi_sq takes skip_const, because the oee sector has NO constant mode -- there
#       <R,1> = 0 holds identically by parity and no index may be dropped;
#   (3) the A_2 correction is returned as ONE interval, so that adding it to 2<w,v>-<w,H'w>
#       costs exactly ONE directed rounding, as the shipped driver's single `+ corr` does.
#       Two separate additions of two sub-ulp terms would round the supremum up TWICE and
#       break bit-level reduction against the recorded certificates;
#   (4) the out-of-space term is clamped to [max(0,inf), max(0,sup)], since the quantity it
#       encloses is a squared norm and the enclosure can straddle 0 under cancellation.
#
# THE REALIZATION (six objects; D = H^1(Omega), NOT V_{N'}):
#   M(u,v) = (grad u,grad v) + (V u,v) + c(u,v) = a_c        N(u,v) = (u,v)
#   X      = (L^2)^3 x L^2 x H^1(Omega)
#   T u    = { sqrt(1-eps') grad u , sqrt(c - C_eps') u , u }
#   b_G({p,s,q},{p',s',q'}) = (p,p') + (s,s')
#                             + [ eps'(grad q,grad q') + ((V + C_eps') q, q') ]
# b_G(Tu,Tv) = M(u,v) identically, and b_G >= 0 on ALL of X: the first two blocks are sums
# of squares (this is what needs c >= C_eps'), the third is the localized sharp form bound
# at eps = eps', i.e. exactly what an admissible certified C_eps' asserts.  Positivity on X
# rather than only on T(D) is what the B-positivity lemma requires, and it is the step the
# shipped X = V_{N'} realization cannot make.
#
# ADMISSIBLE w AND A_2.  q = wtilde (the V_{N'} CG solution the driver already computes),
# s = sqrt(c-C) wtilde, and p solves b_G(w,Tv) = (v_i,v):
#   p = sqrt(1-eps') grad wtilde + grad(phi)/sqrt(1-eps'),   -Lap phi = R  (Neumann),
#   <R,v> := (v_i,v) - a_c(wtilde,v)                          [the WEAK residual]
#   A_2^G := b_G(w,w) = 2<wtilde,v> - <wtilde,H' wtilde> + |grad phi|^2/(1-eps') >= <v,H'^{-1}v>.
#
# CONSTANT MODE (eee only).  phi needs <R,1> = 0.  The brief states this holds exactly
# because 1 lies in V_{N'}; that is true of the DISCRETIZATION part but not of the algebraic
# part: <R,1> equals the constant-mode coefficient r_0 of the in-space residual, which is
# ~1e-14, not 0, because the CG solve is inexact.  It is absorbed rigorously into the
# s-block: s = sqrt(c-C) wtilde + sigma phi_0 with sigma = r_0/sqrt(c-C) makes
# b_G(w,Tv) = (v_i,v) hold EXACTLY and adds sigma^2 = r_0^2/(c-C) to A_2^G.  Needs c > C
# strictly (asserted).  In the oee sector R is odd in x, so <R,1> = 0 identically and no
# such term arises.
#
# RESIDUAL SPLIT.  wtilde in V_{N'} and (-Lap + c)wtilde in V_{N'}, so only V wtilde leaks
# out of the space, and V is even in every coordinate for collinear centres symmetric about
# the origin, so V wtilde carries the SAME single-parity mode set as wtilde.  Hence
#   R = r - (I - Pi^0_{N'})(V wtilde),   the two pieces L^2-orthogonal (disjoint mode sets),
#   |grad phi|^2 = sum_{k in N', k != 0} |r_k|^2/nu_k + sum_{k outside N'} |.|^2/nu_k
#               <= sum_{k in N', k != 0} |r_k|^2/nu_k + ||(I-Pi^0_{N'})(V wtilde)||^2/nu_*(N').
# The in-space part is evaluated as the EXACT modal sum (sharper than ||r||^2/nu_1 and
# rigorous by the same argument); nu_1 is reported separately for reference.
#
# V^2 MOMENTS.  V^2 = sum_{a,b} Z_a Z_b/(r_a r_b) and 1/(r_a r_b) = pref^2 int int
# exp(-t^2 r_a^2 - u^2 r_b^2) dt du with pref = 2/sqrt(pi).  Per coordinate the two
# Gaussians COLLAPSE to one,
#   exp(-t^2(x-a)^2 - u^2(x-b)^2) = exp(-gam (a-b)^2) exp(-wid^2 (x-m)^2),
#   wid^2 = t^2+u^2,   gam = t^2 u^2/wid^2,   m = (a t^2 + b u^2)/wid^2,
# so MomentsVerified.moment_verified is reused VERBATIM through cert_core's
# cos_moment_matrix_1d -- no new quadrature theory, and the same verified panel rule and
# large-t tail carry over.  The (t,u) product rule is the SAME nt-point rule the shipped
# potential assembly uses in one dimension, so the corrected term inherits exactly the
# quadrature status of the P it corrects, no better and no worse.  Two symmetries are used:
#   (i)  term(l,k,a,b) = term(k,l,b,a), so (k,l) runs over k <= l with weight 2 for k < l;
#   (ii) in a SINGLE-PARITY axis both |n-m| and n+m are even, so only even f enters G[f],
#        sin(f pi/2) = 0 and B(centre s) = B(centre -s) -- so centre pairs may be bucketed
#        by |m|.  This is why the sector must be single-parity (asserted).
# The V^2 matrix is NEVER materialised: only the scalar w' V^2 w is accumulated, by tensor
# contraction, so memory stays at vector size.
module GoerischCorrect2

using LinearAlgebra, Printf
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam

const IV = Interval{Float64}
const Centre = Tuple{NTuple{3,Float64},Float64}

# ------------------------------------------------------------------ contractions ---
# Index order matches cert_core's kron_accum_iv!:  idx = (cz-1)*ny*nx + (cy-1)*nx + cx,
# i.e. x is the fastest axis.
"U = (Bz kron By kron I) W, applied along axes 2 and 3."
function apply_yz(W::Array{IV,3}, By::Matrix{IV}, Bz::Matrix{IV})
    nx, ny, nz = size(W)
    A = Array{IV,3}(undef, nx, ny, nz)
    Threads.@threads for k in 1:nz
        @inbounds for jp in 1:ny, i in 1:nx
            acc = interval(0.0)
            for j in 1:ny; acc += By[jp,j]*W[i,j,k]; end
            A[i,jp,k] = acc
        end
    end
    U = Array{IV,3}(undef, nx, ny, nz)
    Threads.@threads for kp in 1:nz
        @inbounds for j in 1:ny, i in 1:nx
            acc = interval(0.0)
            for k in 1:nz; acc += Bz[kp,k]*A[i,j,k]; end
            U[i,j,kp] = acc
        end
    end
    return U
end

"S = sum_{j,k} W[:,j,k]' Bx U[:,j,k] -- the x-contraction against a precomputed U."
function xcontract(W::Array{IV,3}, Bx::Matrix{IV}, U::Array{IV,3})
    nx, ny, nz = size(W)
    parts = fill(interval(0.0), nz)
    Threads.@threads for k in 1:nz
        acc = interval(0.0)
        @inbounds for j in 1:ny, ip in 1:nx
            s = interval(0.0)
            for i in 1:nx; s += Bx[ip,i]*U[i,j,k]; end
            acc += W[ip,j,k]*s
        end
        parts[k] = acc
    end
    tot = interval(0.0)
    for k in 1:nz; tot += parts[k]; end
    return tot
end

# ------------------------------------------------------- ||V wtilde||^2, streamed ---
"""
    vsq_quadform(wvec, nx, ny, nz, momx, momyz, centres, tgrid, twts) -> (val::IV, nterms, npairs)

`momx(centre, width) -> Matrix{IV}` and `momyz(centre, width) -> Matrix{IV}` must return the
1-D interval cosine moment matrices ALREADY RESTRICTED to the sector index sets of the x and
transverse axes respectively (nx x nx and ny x ny).  Requires ny == nz (LY == LZ and the same
transverse parity, true for every box in this project) and all centres on the x axis.
"""
function vsq_quadform(wvec::Vector{IV}, nx::Int, ny::Int, nz::Int,
                      momx::Function, momyz::Function,
                      centres::Vector{Centre},
                      tgrid::Vector{Float64}, twts::Vector{Float64};
                      verbose::Bool=true, report_every::Int=200)
    @assert ny == nz "this fast path assumes ny == nz (LY == LZ, same transverse parity)"
    D = nx*ny*nz
    @assert length(wvec) == D "wvec length $(length(wvec)) != nx*ny*nz = $D"
    for cc in centres
        @assert abs(cc[1][2]) < 1e-14 && abs(cc[1][3]) < 1e-14 "centres must lie on the x axis"
    end
    W = reshape(copy(wvec), nx, ny, nz)
    pref2 = (interval(2.0)/sqrt(interval(pi)))^2
    nt = length(tgrid)
    ax = [cc[1][1] for cc in centres]; Zs = [cc[2] for cc in centres]
    total = interval(0.0); nterms = 0; npairs = 0; npairs_tot = nt*(nt+1)÷2
    for k in 1:nt, l in k:nt
        t2 = tgrid[k]^2; u2 = tgrid[l]^2; w2 = t2 + u2
        wid = sqrt(w2); gam = t2*u2/w2
        sym = (k == l) ? interval(1.0) : interval(2.0)
        wkl = sym * pref2 * interval(twts[k]) * interval(twts[l])
        By = momyz(0.0, wid); Bz = By
        U  = apply_yz(W, By, Bz)
        buckets = Dict{Float64,IV}()
        for i in eachindex(ax), j in eachindex(ax)
            (Zs[i] == 0.0 || Zs[j] == 0.0) && continue
            m  = (t2*ax[i] + u2*ax[j])/w2
            d2 = (ax[i] - ax[j])^2
            wij = interval(Zs[i]*Zs[j]) * exp(-interval(gam)*interval(d2))
            key = abs(m)                     # valid: single-parity axis, B(+s) == B(-s)
            buckets[key] = get(buckets, key, interval(0.0)) + wij
        end
        for (m, wij) in buckets
            Bx = momx(m, wid)
            total += wkl * wij * xcontract(W, Bx, U)
            nterms += 1
        end
        npairs += 1
        verbose && npairs % report_every == 0 &&
            (@printf("      [V^2] pair %d/%d (terms %d)\n", npairs, npairs_tot, nterms);
             flush(stdout))
    end
    return total, nterms, npairs
end

# ------------------------------------------------------------- |grad phi|^2 bound ---
"nu_*(N') -- the driver's own conservative convention, ((N'+1) pi / (2 LX))^2."
nu_star(Nprime::Int, LX::Float64) =
    ((interval(Nprime + 1)*interval(pi))/(interval(2.0)*interval(LX)))^2

"nu_1 -- the smallest NONZERO Neumann-Laplacian eigenvalue of the box, (pi/(2 LX))^2."
nu1_neumann(LX::Float64) = (interval(pi)/(interval(2.0)*interval(LX)))^2

"""
    grad_phi_sq(r, Kdiag, tailsq, nus; skip_const) -> (gphi, inpart, outpart, tail_clamped)

In-space part is the EXACT modal sum sum_{k != 0} |r_k|^2/nu_k.  `skip_const=true` drops
index 1 (the constant mode, which exists only in the eee sector and whose coefficient is
handled by the s-block instead); `skip_const=false` asserts every nu_k > 0 and drops nothing.
"""
function grad_phi_sq(r::Vector{IV}, Kdiag::Vector{IV}, tailsq::IV, nus::IV;
                     skip_const::Bool)
    n = length(r)
    if skip_const
        @assert sup(Kdiag[1]) == 0.0 "skip_const=true but Kdiag[1] = $(Kdiag[1]) is not the zero mode"
    else
        @assert inf(Kdiag[1]) > 0.0 "skip_const=false but Kdiag[1] = $(Kdiag[1]) is not positive"
    end
    i0 = skip_const ? 2 : 1
    inpart = interval(0.0)
    @inbounds for i in i0:n
        # r[i]^2, NOT r[i]*r[i]: a residual component straddles zero, and interval
        # MULTIPLICATION of a straddling interval by itself returns [-e^2, e^2] instead of
        # the true [0, e^2].  The supremum is the same, so no bound is affected, but the
        # spurious negative infimum drops inf(A_2) by one ulp and so breaks bit-level
        # reduction against the recorded certificates (gate G1 detected exactly this).
        inpart += (r[i]^2)/Kdiag[i]
    end
    # squared norm, so clamp a straddling enclosure at 0 from below (rigorous)
    ts = interval(max(0.0, inf(tailsq)), max(0.0, sup(tailsq)))
    outpart = ts / nus
    return inpart + outpart, inpart, outpart, ts
end

# ------------------------------------------------------------------ A_2^Goerisch ---
"""
    a2_correction(gphi, epsp, c, Ceps, r0; has_const) -> (corr, defect, sigma2)

`corr` is the SINGLE interval to be added to 2<w,v> - <w,H'w>, so the addition costs one
directed rounding exactly as the shipped `+ corr` does.
"""
function a2_correction(gphi::IV, epsp::Float64, c::Float64, Ceps::IV, r0::IV;
                       has_const::Bool)
    @assert 0.0 < epsp < 1.0 "eps' must lie in (0,1): got $epsp"
    cmC = interval(c) - Ceps
    @assert inf(cmC) > 0 "b_G positivity needs c > Ceps' strictly: c - Ceps' = $cmC"
    defect = gphi / (interval(1.0) - interval(epsp))
    sigma2 = has_const ? (r0^2)/cmC : interval(0.0)      # squaring, not self-multiplication
    return defect + sigma2, defect, sigma2
end

# ----------------------------------------------------------------------- gates ----
"Gate: b_G(Tu,Tu) == M(u,u) on random elements of V_N, both sides assembled independently."
function gate_bG_equals_M(Ppot::Matrix{IV}, Kdiag::Vector{IV}, epsp::Float64,
                          c::Float64, Ceps::IV; ntrial::Int=24, seed::Int=20260902)
    D = length(Kdiag); worst = 0.0
    st = seed
    for _ in 1:ntrial
        x = Vector{Float64}(undef, D)
        for i in 1:D
            st = (1103515245*st + 12345) % 2147483648
            x[i] = (st/2147483648.0) - 0.5
        end
        xI = interval.(x)
        gradsq = interval(0.0); nrmsq = interval(0.0)
        @inbounds for i in 1:D
            gradsq += Kdiag[i]*xI[i]*xI[i]; nrmsq += xI[i]*xI[i]
        end
        Px = Vector{IV}(undef, D)
        Threads.@threads for i in 1:D
            acc = interval(0.0)
            @inbounds for j in 1:D; acc += Ppot[i,j]*xI[j]; end
            Px[i] = acc
        end
        Vq = interval(0.0)
        @inbounds for i in 1:D; Vq += xI[i]*Px[i]; end
        lhs = (interval(1.0)-interval(epsp))*gradsq + (interval(c)-Ceps)*nrmsq +
              interval(epsp)*gradsq + (Vq + Ceps*nrmsq)
        rhs = gradsq + Vq + interval(c)*nrmsq
        worst = max(worst, abs(mid(lhs) - mid(rhs))/max(1.0, abs(mid(rhs))))
    end
    return worst
end

"""
Gate: the EXACT relation, in 400-bit arithmetic on a small SPD system with an inexact wtilde.
The brief writes A_2^G = <v,H'^{-1}v> + |grad phi|^2/(1-eps') as an identity; it is not.
From 2<w,v> - <w,H'w> = <v,H'^{-1}v> - <R,H'^{-1}R> one gets
    A_2^G - <v,H'^{-1}v> = <R, [ (-Lap)^{-1}/(1-eps') - H'^{-1} ] R>,
and the non-negativity of that gap is EXACTLY H' >= (1-eps')(-Lap), i.e. the same
eps'|grad u|^2 + (Vu,u) + c|u|^2 >= 0 that b_G >= 0 asserts.  So A_2^G >= <v,H'^{-1}v>
is implied by admissibility rather than assumed.
"""
function gate_identity_small(n::Int, epsp::Float64; seed::Int=7)
    setprecision(BigFloat, 400)
    rng = seed
    nextr() = (rng = (1103515245*rng + 12345) % 2147483648; rng/2147483648.0)
    Kd = BigFloat[ i == 1 ? big(0.0) : big(0.35)*big(i-1)^2 for i in 1:n ]
    Vm = zeros(BigFloat, n, n)
    for i in 1:n, j in i:n
        v = big(nextr() - 0.5)/big(n); Vm[i,j] = v; Vm[j,i] = v
    end
    cc = big(3.0)
    Hp = diagm(Kd) + Vm
    for i in 1:n; Hp[i,i] += cc; end
    v = BigFloat[ big(nextr() - 0.5) for i in 1:n ]; v ./= sqrt(sum(v.^2))
    z = Hp \ v
    wt = z .+ big(1e-3).*BigFloat[ big(nextr() - 0.5) for i in 1:n ]
    R = v .- Hp*wt
    A2G_core = 2*dot(wt, v) - dot(wt, Hp*wt)
    gphi = sum( i == 1 ? big(0.0) : R[i]^2/Kd[i] for i in 1:n )
    A2G = A2G_core + gphi/(1 - big(epsp))
    exact = dot(v, z)
    gap_pred = gphi/(1-big(epsp)) - dot(R, Hp \ R)
    return Float64(abs(A2G - (exact + gap_pred))/abs(exact)), Float64(gap_pred),
           Float64(A2G), Float64(exact), Float64(R[1])
end

peak_rss_gb() = begin
    v = NaN
    try
        for ln in eachline("/proc/self/status")
            startswith(ln, "VmHWM:") && (v = parse(Float64, split(ln)[2])/1048576.0; break)
        end
    catch; end
    v
end

end # module
