# cert_core.jl -- box- and system-parameterised CERTIFIED (interval) chain.
# =============================================================================
# The rigorous counterpart of the float sweep.  This is a PARAMETERISATION of the
# bundle's validated Stage-2 modules, not a reimplementation: every rigorous step
# (moment enclosure, Lehmann-Behnke Galerkin enclosure, Goerisch A2 with its
# residual correction, the Moebius transform, the Dirichlet Rayleigh quotient) is
# the bundle's own algebra in the bundle's own operation order.  moments_verified.jl
# is used byte-for-byte as shipped (the omega1 and omega2 copies differ only in a
# default npanel, which is always passed explicitly here).
#
# Deliberate differences from the published scripts, all disclosed:
#  (D1) BOX AND CENTRES ARE ARGUMENTS.  The published modules hardcode
#       `const LX,LY,LZ` and `const AX = 2.0` with exactly two centres.  Here the
#       box and an arbitrary centre list are passed, using the grouping algebra of
#       the validated assembly_multicentre.jl (centres sharing (ay,az,Z) have their
#       Bx factors summed BEFORE the Kronecker product).  For the two H2+ centres
#       and a published box this reduces to the published code bit-identically.
#  (D2) npanel AND t_star ARE ARGUMENTS, set by the rules established in the float
#       sweep: npanel = round(4.8*LX) (which reproduces the published defaults, 48
#       at LX=10 and 96 at LX=20), and t_star raised from 1.0 in 0.25 steps until
#       the large-t moment tail bound is <= 1e-13.
#  (D3) THE INTERVAL KRONECKER ACCUMULATION IS FUSED AND THREADED.  The published
#       code evaluates `P .+= wk .* kron(Bz, kron(By, Bx))`, which allocates a
#       second D x D interval temporary (76 GB at D=68921) and runs on one thread.
#       Here the same expression is evaluated entry by entry in the same
#       association order,
#           P[i,j] <- P[i,j] + wkZ * (Bz[iz,jz] * (By[iy,jy] * Bx[ix,jx])),
#       in parallel over columns, one writer per entry.  cert_selftest.jl VERIFIES
#       bitwise on both endpoints that threaded == serial == the published
#       kron form, and the production driver REFUSES to run if it does not.
#       This is what makes the runs affordable: published timings scale as
#       5.87e-6 * D^2 s single-threaded, i.e. ~7.8 h for one N'=80 assembly.
#  (D4) form_H IS IN PLACE.  The published `form_H` does `H = copy(P)`, a second
#       D x D interval matrix.  Here the diagonal is modified in place after
#       SAVING the original diagonal entries, and restored by assignment (never
#       by interval subtraction, which would widen).
#
# NOT changed: the approximate eigensolver stays `Arpack.eigs(Symmetric(Hm); nev,
# which=:SR)` exactly as the published scripts use it.  The provider doc's Arpack
# warning concerns the `sigma=` shift-invert form, which is not used here; only the
# subsequent Lehmann-Behnke step is rigorous, so this choice affects bound quality
# and never bound validity.
# =============================================================================
module CertCore

using LinearAlgebra, SparseArrays, Printf
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
include("moments_verified.jl")
using .MomentsVerified
using Arpack, IterativeSolvers
import Veigs: lehmann_behnke

const IV = Interval{Float64}
const Centre = Tuple{NTuple{3,Float64},Float64}

# Set to false (env CERT_THREADED=0) to run every interval accumulation serially.
# Only cert_selftest.jl passing authorises the threaded path.
const THREADED = Ref(true)

h2plus_centres(ax::Float64=2.0) = Centre[((-ax,0.0,0.0),1.0), ((ax,0.0,0.0),1.0)]
h3plus_centres(dp::Float64)     = Centre[((-dp,0.0,0.0),1.0), ((0.0,0.0,0.0),1.0),
                                         ((dp,0.0,0.0),1.0)]

# ------------------------------------------------------ rules from the sweep --
npanel_rule(LX::Float64) = round(Int, 4.8 * LX)

function moment_tail(L::Float64, s::Float64, t::Float64)
    dL = L + s; dR = L - s
    (dL <= 0 || dR <= 0) && return Inf
    return exp(-t*t*dL*dL)/(t*t*dL) + exp(-t*t*dR*dR)/(t*t*dR)
end
cell_tail(cen::Vector{Centre}, LX, LY, LZ, ts::Float64) =
    maximum(max(moment_tail(LX, c[1][1], ts), moment_tail(LY, c[1][2], ts),
                moment_tail(LZ, c[1][3], ts)) for c in cen if c[2] != 0.0)
function choose_t_star(cen::Vector{Centre}, LX, LY, LZ; tol::Float64=1e-13)
    ts = 1.0
    while ts < 8.0 && cell_tail(cen, LX, LY, LZ, ts) > tol; ts += 0.25; end
    return ts, cell_tail(cen, LX, LY, LZ, ts)
end

# sigma(Omega) = inf_{x notin Omega} V, attained at a face centre for centres on
# the x-axis.  Reproduces the recorded H2+ values -0.242535625 / -0.124034735 and
# all 16 recorded H3^2+ confinement rows to 4.8e-10.
function sigma_conf(cen::Vector{Centre}, LX, LY, LZ)
    V(x,y,z) = -sum(c[2]/sqrt((x-c[1][1])^2 + (y-c[1][2])^2 + (z-c[1][3])^2) for c in cen)
    return min(V(LX,0.0,0.0), V(-LX,0.0,0.0), V(0.0,LY,0.0), V(0.0,-LY,0.0),
               V(0.0,0.0,LZ), V(0.0,0.0,-LZ))
end

# ------------------------------------------------------------- basis pieces --
sector_idx(N::Int, p::Symbol) = p === :even ? collect(1:2:(N+1)) : collect(2:2:(N+1))
odd_modes(mmax::Int) = collect(1:2:mmax)

function build_t_grid(nt::Int)
    k = 1:(nt-1); beta = collect(k ./ sqrt.(4.0 .* k.^2 .- 1.0))
    J = SymTridiagonal(zeros(nt), beta); vals, vecs = eigen(J)
    x = vals; w = 2.0 .* (vecs[1, :] .^ 2)
    u = 0.5 .* (x .+ 1); wu = 0.5 .* w
    return u ./ (1 .- u), wu ./ (1 .- u).^2
end

# cosine (Neumann) interval moment matrix -- identical algebra to assembly_verified.jl
function cos_moment_matrix_1d(N::Int, L::Float64, s::Float64, t::Float64,
                              nodes::Vector{IV}, wts::Vector{IV};
                              npanel::Int, t_star::Float64)
    n1 = N+1; fmax = 2N
    G = Vector{IV}(undef, fmax+1)
    for f in 0:fmax
        G[f+1] = MomentsVerified.moment_verified(f*pi/(2L), t, L, s;
                     t_star=t_star, nodes=nodes, wts=wts, npanel=npanel)
    end
    Nn = [ (n==0) ? interval(1.0)/sqrt(interval(2.0)*interval(L)) :
                    interval(1.0)/sqrt(interval(L)) for n in 0:N ]
    B = Matrix{IV}(undef, n1, n1); half = interval(0.5)
    @inbounds for n in 0:N, m in 0:N
        B[n+1,m+1] = half * Nn[n+1] * Nn[m+1] * (G[abs(n-m)+1] + G[n+m+1])
    end
    return B
end

# sine (Dirichlet) interval moment matrix -- identical algebra to dirichlet_verified.jl
function sin_moment_matrix_1d(mmax::Int, L::Float64, s::Float64, t::Float64,
                              nodes::Vector{IV}, wts::Vector{IV};
                              npanel::Int, t_star::Float64)
    fmax = 2*mmax
    G = Vector{IV}(undef, fmax+1)
    for f in 0:fmax
        G[f+1] = MomentsVerified.moment_verified(f*pi/(2L), t, L, s;
                     t_star=t_star, nodes=nodes, wts=wts, npanel=npanel)
    end
    Ninv = interval(1.0)/interval(L); half = interval(0.5)
    B = Matrix{IV}(undef, mmax, mmax)
    @inbounds for n in 1:mmax, m in 1:mmax
        B[n,m] = half*Ninv*(G[abs(n-m)+1] - G[n+m+1])
    end
    return B
end

function group_centres(centres::Vector{Centre})
    ks = Tuple{Float64,Float64,Float64}[]; lists = Vector{Vector{Float64}}()
    for ((ax,ay,az), Z) in centres
        Z == 0.0 && continue
        k = (ay, az, Z); i = findfirst(==(k), ks)
        if i === nothing; push!(ks, k); push!(lists, Float64[ax])
        else; push!(lists[i], ax); end
    end
    return [(ks[i][1], ks[i][2], ks[i][3], lists[i]) for i in eachindex(ks)]
end

# ------------------------- (D3) fused, threaded interval kron accumulation ----
function kron_accum_iv!(P::Matrix{IV}, wkZ::IV,
                        Bx::Matrix{IV}, By::Matrix{IV}, Bz::Matrix{IV};
                        threaded::Bool=THREADED[])
    nx = size(Bx,1); ny = size(By,1); nz = size(Bz,1)
    ncol = nx*ny*nz
    function body(jcol)
        j0 = jcol - 1
        jx = (j0 % nx) + 1
        jy = ((j0 ÷ nx) % ny) + 1
        jz = (j0 ÷ (nx*ny)) + 1
        @inbounds for iz in 1:nz
            bz = Bz[iz,jz]
            for iy in 1:ny
                byy = By[iy,jy]
                base = ((iz-1)*ny + (iy-1))*nx
                for ix in 1:nx
                    P[base+ix, jcol] = P[base+ix, jcol] + wkZ*(bz*(byy*Bx[ix,jx]))
                end
            end
        end
    end
    if threaded
        Threads.@threads for jcol in 1:ncol; body(jcol); end
    else
        for jcol in 1:ncol; body(jcol); end
    end
    return P
end

# ------------------------------------------------ Neumann interval assembly --
function assemble_PK(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                     centres::Vector{Centre}, LX::Float64, LY::Float64, LZ::Float64,
                     nt::Int=48, npanel::Int, t_star::Float64,
                     nodes::Vector{IV}, wts::Vector{IV},
                     threaded::Bool=THREADED[], verbose::Bool=true)
    t, wt = build_t_grid(nt)
    Ix, Iy, Iz = sector_idx(N,px), sector_idx(N,py), sector_idx(N,pz)
    nx, ny, nz = length(Ix), length(Iy), length(Iz)
    D = nx*ny*nz
    pref = interval(2.0)/sqrt(interval(pi))
    groups = group_centres(centres)
    P = fill(interval(0.0), D, D)
    for k in 1:nt
        wk = -pref * interval(wt[k])
        for (ay, az, Z, axs) in groups
            Bx = cos_moment_matrix_1d(N, LX, axs[1], t[k], nodes, wts; npanel=npanel, t_star=t_star)
            for jj in 2:length(axs)
                Bx = Bx .+ cos_moment_matrix_1d(N, LX, axs[jj], t[k], nodes, wts; npanel=npanel, t_star=t_star)
            end
            Bxs = Bx[Ix,Ix]
            By = cos_moment_matrix_1d(N, LY, ay, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iy,Iy]
            Bz = cos_moment_matrix_1d(N, LZ, az, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iz,Iz]
            wkZ = Z == 1.0 ? wk : wk * interval(Z)
            kron_accum_iv!(P, wkZ, Bxs, By, Bz; threaded=threaded)
        end
        verbose && k % 8 == 0 && (@printf("      [iv-asm] t %d/%d\n", k, nt); flush(stdout))
    end
    nux = (interval(pi)/(interval(2.0)*interval(LX)))^2
    nuy = (interval(pi)/(interval(2.0)*interval(LY)))^2
    nuz = (interval(pi)/(interval(2.0)*interval(LZ)))^2
    mx = Ix .- 1; my = Iy .- 1; mz = Iz .- 1
    Kdiag = Vector{IV}(undef, D)
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        idx = (cz-1)*ny*nx + (cy-1)*nx + cx
        Kdiag[idx] = nux*interval(mx[cx]^2) + nuy*interval(my[cy]^2) + nuz*interval(mz[cz]^2)
    end
    return P, Kdiag, D
end

# ---------------------------------------------- Dirichlet interval assembly --
function assemble_dirichlet(N::Int; centres::Vector{Centre},
                            LX::Float64, LY::Float64, LZ::Float64,
                            nt::Int=48, npanel::Int, t_star::Float64,
                            nodes::Vector{IV}, wts::Vector{IV},
                            threaded::Bool=THREADED[], verbose::Bool=true)
    mmax = N+1
    t, wt = build_t_grid(nt)
    Ix = odd_modes(mmax); Iy = Ix; Iz = Ix
    nx = length(Ix); D = nx^3
    pref = interval(2.0)/sqrt(interval(pi))
    groups = group_centres(centres)
    P = fill(interval(0.0), D, D)
    for k in 1:nt
        wk = -pref * interval(wt[k])
        for (ay, az, Z, axs) in groups
            Bx = sin_moment_matrix_1d(mmax, LX, axs[1], t[k], nodes, wts; npanel=npanel, t_star=t_star)
            for jj in 2:length(axs)
                Bx = Bx .+ sin_moment_matrix_1d(mmax, LX, axs[jj], t[k], nodes, wts; npanel=npanel, t_star=t_star)
            end
            Bxs = Bx[Ix,Ix]
            By = sin_moment_matrix_1d(mmax, LY, ay, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iy,Iy]
            Bz = sin_moment_matrix_1d(mmax, LZ, az, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iz,Iz]
            wkZ = Z == 1.0 ? wk : wk * interval(Z)
            kron_accum_iv!(P, wkZ, Bxs, By, Bz; threaded=threaded)
        end
        verbose && k % 8 == 0 && (@printf("      [iv-dir] t %d/%d\n", k, nt); flush(stdout))
    end
    nux = (interval(pi)/(interval(2.0)*interval(LX)))^2
    nuy = (interval(pi)/(interval(2.0)*interval(LY)))^2
    nuz = (interval(pi)/(interval(2.0)*interval(LZ)))^2
    @inbounds for cz in 1:nx, cy in 1:nx, cx in 1:nx
        idx = (cz-1)*nx*nx + (cy-1)*nx + cx
        P[idx,idx] = P[idx,idx] + nux*interval(Ix[cx]^2) + nuy*interval(Iy[cy]^2) + nuz*interval(Iz[cz]^2)
    end
    return P, D, nx
end

# ----------------------------------- (D4) in-place H = P + ks*K + sigma I ----
function form_H_inplace!(P::Matrix{IV}, Kdiag::Vector{IV};
                         kin_scale=interval(1.0), sigma=interval(0.0))
    D = length(Kdiag)
    saved = Vector{IV}(undef, D)
    ks = kin_scale isa IV ? kin_scale : interval(kin_scale)
    sg = sigma isa IV ? sigma : interval(sigma)
    @inbounds for i in 1:D
        saved[i] = P[i,i]
        P[i,i] = P[i,i] + ks*Kdiag[i] + sg
    end
    return saved
end
restore_diag!(P::Matrix{IV}, saved::Vector{IV}) =
    (@inbounds for i in eachindex(saved); P[i,i] = saved[i]; end; P)

# --------------------------------------------------- threaded interval ops ---
# Row i sums over j in increasing order, exactly as the published imatvec.
function imatvec(M::Matrix{IV}, x::Vector{IV}; threaded::Bool=THREADED[])
    D = length(x); y = Vector{IV}(undef, D)
    function body(i)
        acc = interval(0.0)
        @inbounds for j in 1:D; acc += M[i,j]*x[j]; end
        y[i] = acc
    end
    if threaded; Threads.@threads for i in 1:D; body(i); end
    else; for i in 1:D; body(i); end; end
    return y
end
idot(a::Vector{IV}, b::Vector{IV}) =
    (s = interval(0.0); @inbounds for i in eachindex(a); s += a[i]*b[i]; end; s)

function mid_matrix(M::Matrix{IV})
    D = size(M,1); A = Matrix{Float64}(undef, D, D)
    Threads.@threads for j in 1:D
        @inbounds for i in 1:D; A[i,j] = mid(M[i,j]); end
    end
    return A
end

# ------------------------------------------- verified Galerkin enclosures ----
# Verbatim from Pipeline2b.lb_ground / LGOEE.lb_index.
function lb_index(Hm::Matrix{Float64}, Hi::Matrix{IV}, D::Int, k::Int,
                  rho::Float64, sigma::Float64; nev::Int=4)
    vals, vecs = eigs(Symmetric(Hm); nev=nev, which=:SR, maxiter=3000)
    p = sortperm(real(vals)); vals = real(vals[p]); vecs = real(vecs[:,p])
    V = reshape(vecs[:,k], D, 1)
    Bi = interval.(sparse(1.0I, D, D))
    lg = lehmann_behnke(Hi, Bi, vals, V, rho, sigma, 1.0, k, k; do_shift=true)
    return lg[1], vals, vecs
end

# --------------------------------------------------- sharp-L (Stage-A) -------
nu_star(N::Int, LX::Float64) = (interval(N+1)*interval(pi)/(interval(2.0)*interval(LX)))^2
function sharp_L_iv(m1::IV, mk::IV, N::Int, eps::Float64, Ceps::IV, LX::Float64)
    e = interval(eps); one_ = interval(1.0); sg = Ceps
    ms1 = m1 + sg
    inf(ms1) > 0 || return nothing
    etaV = sqrt(e/(one_-e) + Ceps/ms1)
    g = etaV/sqrt(ms1)
    Ch2 = (one_ + g*g)/((one_-e)*nu_star(N, LX))
    msk = mk + sg
    return msk/(one_ + Ch2*msk) - sg
end

# ------------------------- (D7) midpoint-on-the-fly operator + own CG --------
# The published LG scripts materialise `Hm = mid.(Hc)`, a full D x D FLOAT matrix,
# purely to run CG and to slice the primary block out of.  At D=68921 that is a
# 38 GB allocation on top of the 76 GB interval operator.  /dev/shm on this host is
# currently holding ~85 GB of another user's tmpfs, leaving ~159 GB, and the
# Lehmann-Behnke path already needs several dense interval copies -- so the float
# midpoint is the allocation to remove.  MidOp computes mid(P[i,j]) during the
# matvec instead (threaded over rows, one streaming pass over P per iteration), and
# cg_mid is a plain CG against it.
#
# This is a FLOAT preconditioning step: w only has to be a good approximate solve.
# Rigour comes from the INTERVAL residual r = v - H'w computed afterwards, so a w
# from a different-but-converged CG changes the bound in its last digits and can
# never invalidate it.
struct MidOp
    P::Matrix{IV}
end
Base.size(A::MidOp) = size(A.P)
Base.size(A::MidOp, i::Int) = size(A.P, i)
Base.eltype(::MidOp) = Float64

function midmul!(y::Vector{Float64}, A::MidOp, x::Vector{Float64})
    P = A.P; D = length(x)
    Threads.@threads for i in 1:D
        s = 0.0
        @inbounds for j in 1:D; s += mid(P[i,j])*x[j]; end
        y[i] = s
    end
    return y
end

"""
    cg_mid(A::MidOp, b; reltol=1e-12, maxiter=8000) -> (x, iters, relres)

Conjugate gradients for the SPD system mid(P) x = b, matching the published
`cg!(w, Symmetric(Hm), v; reltol=1e-12, maxiter=8000)` call but without
materialising the float midpoint matrix.
"""
function cg_mid(A::MidOp, b::Vector{Float64}; reltol::Float64=1e-12, maxiter::Int=8000)
    n = length(b)
    x = zeros(Float64, n); r = copy(b); p = copy(r); Ap = similar(r)
    rs = dot(r, r); nb = sqrt(rs); it = 0
    for i in 1:maxiter
        it = i
        midmul!(Ap, A, p)
        denom = dot(p, Ap)
        denom <= 0 && break
        alpha = rs/denom
        axpy!(alpha, p, x); axpy!(-alpha, Ap, r)
        rs_new = dot(r, r)
        sqrt(rs_new) <= reltol*nb && break
        p .= r .+ (rs_new/rs) .* p
        rs = rs_new
    end
    midmul!(Ap, A, x)
    return x, it, norm(Ap .- b)/nb
end

# mid of a symmetric sub-block, without materialising the full float midpoint
function mid_subblock(P::Matrix{IV}, idx::Vector{Int})
    m = length(idx); S = Matrix{Float64}(undef, m, m)
    Threads.@threads for jj in 1:m
        j = idx[jj]
        @inbounds for ii in 1:m; S[ii,jj] = mid(P[idx[ii], j]); end
    end
    return S
end

function peak_rss_gb()
    try
        for ln in eachline("/proc/self/status")
            startswith(ln, "VmHWM:") && return parse(Float64, split(ln)[2]) / 1048576.0
        end
    catch; end
    return NaN
end

end # module
