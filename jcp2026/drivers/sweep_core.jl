# sweep_core.jl -- multi-centre float assembly + shift-invert ground solver for
# the (L,N) enclosure-width sweep.
#
# Extends float_core.jl with
#   * MULTI-CENTRE Coulomb assembly (Neumann cosine and Dirichlet sine sectors),
#     using the grouping algebra of assembly_multicentre.jl: centres sharing a
#     transverse position (ay,az) and charge Z have their Bx factors summed
#     BEFORE the Kronecker product, so a collinear arrangement costs exactly one
#     kron per t node regardless of how many nuclei there are.  For
#     centres = [((-2,0,0),1),((2,0,0),1)] this reduces to the two-centre code.
#   * the large-t moment TAIL BOUND, and the t_star selection rule that keeps it
#     below a stated tolerance (the guard against a nucleus near the boundary).
#   * a Cholesky-based SHIFT-INVERT LANCZOS ground solver, replacing
#     eigen(Symmetric(A),1:1): LAPACK dsyevr costs 446 s at D=35937 on this host
#     because the tridiagonal reduction threads poorly, while one dense Cholesky
#     plus ~20 triangular-solve Lanczos steps is several times cheaper.
#     Arpack.eigs(...; sigma=...) is NOT used -- it silently returns the largest
#     eigenvalues on this host.
#
# Still zero external packages: LinearAlgebra + Printf only.
module SweepCore

using LinearAlgebra, Printf, Random
include("float_core.jl")
using .FloatCore
const FC = FloatCore

const Centre = Tuple{NTuple{3,Float64},Float64}

h2plus_centres(ax::Float64=2.0) = Centre[((-ax,0.0,0.0),1.0), ((ax,0.0,0.0),1.0)]
h3plus_centres(dp::Float64)     = Centre[((-dp,0.0,0.0),1.0), ((0.0,0.0,0.0),1.0),
                                         ((dp,0.0,0.0),1.0)]

# ----------------------------------------------------- centre grouping -------
# (ay, az, Z, [ax...]) in first-appearance order; Z == 0 centres dropped.
function group_centres(centres::Vector{Centre})
    keys = Tuple{Float64,Float64,Float64}[]
    lists = Vector{Vector{Float64}}()
    for ((ax,ay,az), Z) in centres
        Z == 0.0 && continue
        k = (ay, az, Z)
        i = findfirst(==(k), keys)
        if i === nothing
            push!(keys, k); push!(lists, Float64[ax])
        else
            push!(lists[i], ax)
        end
    end
    return [(keys[i][1], keys[i][2], keys[i][3], lists[i]) for i in eachindex(keys)]
end

# ------------------------------------------------ large-t moment tail --------
# moment_larget returns the INFINITE-domain closed form; the neglected finite-box
# remainder is bounded by the same two-sided Gaussian tail the interval code adds
# as its uncertainty:
#     tail(L,s,t) = 1/(t^2 dL) exp(-t^2 dL^2) + 1/(t^2 dR) exp(-t^2 dR^2),
#     dL = L+s, dR = L-s.
# It is monotone decreasing in t, so the worst case over the large-t branch is at
# t = t_star.  This is a genuine float error, not just an interval width: it is
# the part of the integral over |x|>L that the closed form wrongly includes.
function moment_tail(L::Float64, s::Float64, t::Float64)
    dL = L + s; dR = L - s
    (dL <= 0 || dR <= 0) && return Inf
    return exp(-t*t*dL*dL)/(t*t*dL) + exp(-t*t*dR*dR)/(t*t*dR)
end

# worst tail over all axes/centres of a cell, at a given t_star
function cell_tail(centres::Vector{Centre}, LX, LY, LZ, t_star::Float64)
    w = 0.0
    for ((ax,ay,az), Z) in centres
        Z == 0.0 && continue
        w = max(w, moment_tail(LX, ax, t_star),
                   moment_tail(LY, ay, t_star), moment_tail(LZ, az, t_star))
    end
    return w
end

# t_star rule: keep the large-t tail at or below `tol`.  t_star = 1.0 (the value
# used for every published H2+ number) whenever that already suffices; otherwise
# raise it in steps of 0.25.  Composite Gauss-Legendre takes over below t_star
# and is exact on the finite box, so raising t_star REMOVES the error rather
# than trading it -- the only cost is more quadrature panels being exercised,
# and with panel half-length h ~ 0.21 and 24 nodes per panel a Gaussian of width
# 1/t_star >= 0.25 is still heavily oversampled.
function choose_t_star(centres::Vector{Centre}, LX, LY, LZ; tol::Float64=1e-13)
    ts = 1.0
    while ts < 8.0 && cell_tail(centres, LX, LY, LZ, ts) > tol
        ts += 0.25
    end
    return ts, cell_tail(centres, LX, LY, LZ, ts)
end

# npanel rule: constant panel half-length in x.  round(4.8*LX) reproduces the
# published defaults exactly -- 48 at LX=10 (Omega1) and 96 at LX=20 (Omega2).
npanel_rule(LX::Float64) = round(Int, 4.8 * LX)

# ------------------------------------------------- multicentre assembly ------
"""
    assemble_neumann_mc(N, px, py, pz; centres, LX, LY, LZ, ...) -> (H, D)

Float Neumann sector Hamiltonian H = P + kin_scale*diag(K) + sigma*I in the
tensor cosine basis.
"""
function assemble_neumann_mc(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                             centres::Vector{Centre},
                             LX::Float64, LY::Float64, LZ::Float64,
                             nt::Int=48, npanel::Int, t_star::Float64=1.0,
                             nodes::Vector{Float64}, wts::Vector{Float64},
                             kin_scale::Float64=1.0, sigma::Float64=0.0,
                             verbose::Bool=false)
    t, wt = FC.build_t_grid(nt)
    Ix, Iy, Iz = FC.sector_idx(N,px), FC.sector_idx(N,py), FC.sector_idx(N,pz)
    nx, ny, nz = length(Ix), length(Iy), length(Iz)
    D = nx*ny*nz
    pref = 2.0/sqrt(pi)
    groups = group_centres(centres)
    H = zeros(Float64, D, D)
    for k in 1:nt
        wk = -pref * wt[k]
        for (ay, az, Z, axs) in groups
            Bx = FC.cos_moment_matrix_1d(N, LX, axs[1], t[k], nodes, wts;
                                         npanel=npanel, t_star=t_star)
            for j in 2:length(axs)
                Bx = Bx .+ FC.cos_moment_matrix_1d(N, LX, axs[j], t[k], nodes, wts;
                                                   npanel=npanel, t_star=t_star)
            end
            Bx = Bx[Ix,Ix]
            By = FC.cos_moment_matrix_1d(N, LY, ay, t[k], nodes, wts;
                                         npanel=npanel, t_star=t_star)[Iy,Iy]
            Bz = FC.cos_moment_matrix_1d(N, LZ, az, t[k], nodes, wts;
                                         npanel=npanel, t_star=t_star)[Iz,Iz]
            FC.kron_accum!(H, wk*Z, Bx, By, Bz)
        end
        verbose && k % 16 == 0 && (@printf("      [asm] %d/%d\n", k, nt); flush(stdout))
    end
    nux = (pi/(2.0*LX))^2; nuy = (pi/(2.0*LY))^2; nuz = (pi/(2.0*LZ))^2
    mx = Ix .- 1; my = Iy .- 1; mz = Iz .- 1
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        idx = ((cz-1)*ny + (cy-1))*nx + cx
        H[idx,idx] += kin_scale*(nux*mx[cx]^2 + nuy*my[cy]^2 + nuz*mz[cz]^2) + sigma
    end
    return H, D
end

"""
    assemble_dirichlet_mc(N; centres, LX, LY, LZ, ...) -> (H, D, nx)

Float Dirichlet Hamiltonian on the even-in-x-y-z sector of the sine basis
psi_n = (1/sqrt(L)) sin(n pi (x+L)/(2L)): the ODD sine modes 1,3,5,... among
1..N+1, giving the same (N/2+1)^3 dimension as the Neumann eee sector.
"""
function assemble_dirichlet_mc(N::Int; centres::Vector{Centre},
                               LX::Float64, LY::Float64, LZ::Float64,
                               nt::Int=48, npanel::Int, t_star::Float64=1.0,
                               nodes::Vector{Float64}, wts::Vector{Float64},
                               verbose::Bool=false)
    mmax = N + 1
    t, wt = FC.build_t_grid(nt)
    Ix = FC.odd_modes(mmax); Iy = Ix; Iz = Ix
    nx = length(Ix); D = nx^3
    pref = 2.0/sqrt(pi)
    groups = group_centres(centres)
    H = zeros(Float64, D, D)
    for k in 1:nt
        wk = -pref * wt[k]
        for (ay, az, Z, axs) in groups
            Bx = FC.sin_moment_matrix_1d(mmax, LX, axs[1], t[k], nodes, wts;
                                         npanel=npanel, t_star=t_star)
            for j in 2:length(axs)
                Bx = Bx .+ FC.sin_moment_matrix_1d(mmax, LX, axs[j], t[k], nodes, wts;
                                                   npanel=npanel, t_star=t_star)
            end
            Bx = Bx[Ix,Ix]
            By = FC.sin_moment_matrix_1d(mmax, LY, ay, t[k], nodes, wts;
                                         npanel=npanel, t_star=t_star)[Iy,Iy]
            Bz = FC.sin_moment_matrix_1d(mmax, LZ, az, t[k], nodes, wts;
                                         npanel=npanel, t_star=t_star)[Iz,Iz]
            FC.kron_accum!(H, wk*Z, Bx, By, Bz)
        end
        verbose && k % 16 == 0 && (@printf("      [asm] %d/%d\n", k, nt); flush(stdout))
    end
    nux = (pi/(2.0*LX))^2; nuy = (pi/(2.0*LY))^2; nuz = (pi/(2.0*LZ))^2
    @inbounds for cz in 1:nx, cy in 1:nx, cx in 1:nx
        idx = ((cz-1)*nx + (cy-1))*nx + cx
        H[idx,idx] += nux*Ix[cx]^2 + nuy*Iy[cy]^2 + nuz*Iz[cz]^2
    end
    return H, D, nx
end

# --------------------------------------- shift-invert Lanczos ground solver --
"""
    ground_si(H, shift; nev=1, mmax=90, tol=1e-13) -> (vals, vecs, iters, info)

Lowest `nev` eigenvalues (and their eigenvectors) of the dense symmetric H.

Method: H - shift*I must be SPD (shift strictly below lambda_1).  Factor it once
with a dense Cholesky IN PLACE (H is destroyed), then run Lanczos with full
reorthogonalisation on B = (H - shift*I)^{-1}, whose LARGEST eigenvalues are the
ones nearest shift from above.  Map back with lambda = shift + 1/theta.

Throws `PosDefException` if the shift is not low enough -- the caller must pick
a safe shift because H is consumed by the factorisation.
"""
function ground_si(H::Matrix{Float64}, shift::Float64;
                   nev::Int=1, mmax::Int=90, tol::Float64=1e-13, seed::Int=20260829)
    D = size(H,1)
    @inbounds for i in 1:D; H[i,i] -= shift; end
    F = cholesky!(Symmetric(H))          # in place; H now holds the factor

    rng = MersenneTwister(seed)
    V = Matrix{Float64}(undef, D, mmax)
    alpha = Float64[]; beta = Float64[]
    q = randn(rng, D); q ./= norm(q)
    V[:,1] .= q
    w = Vector{Float64}(undef, D)
    m = 0; converged = false; theta = Float64[]; Y = zeros(0,0)
    for j in 1:mmax
        m = j
        w .= V[:,j]
        ldiv!(F, w)                       # w = B * v_j
        if j > 1
            axpy!(-beta[j-1], view(V,:,j-1), w)
        end
        a = dot(view(V,:,j), w)
        push!(alpha, a)
        axpy!(-a, view(V,:,j), w)
        # full reorthogonalisation (twice -- cheap at these m, kills drift)
        for _ in 1:2, i in 1:j
            axpy!(-dot(view(V,:,i), w), view(V,:,i), w)
        end
        b = norm(w)
        # Ritz values of the current tridiagonal
        T = SymTridiagonal(copy(alpha), copy(beta))
        es = eigen(T)
        theta = es.values; Y = es.vectors
        if j >= max(nev+4, 12)
            # residual bound for the nev largest theta (= nev lowest lambda)
            ok = true
            for r in 0:(nev-1)
                idx = length(theta) - r
                res = b * abs(Y[end, idx])
                lam = shift + 1.0/theta[idx]
                ok &= (res / (theta[idx]^2) < tol * max(1.0, abs(lam)))
            end
            if ok; converged = true; end
        end
        (converged || b <= 1e-14 || j == mmax) && break
        push!(beta, b)
        V[:,j+1] .= w ./ b
    end
    vals = Float64[]; vecs = Matrix{Float64}(undef, D, nev)
    for r in 0:(nev-1)
        idx = length(theta) - r
        push!(vals, shift + 1.0/theta[idx])
        y = Y[:, idx]
        v = view(V,:,1:m) * y
        v ./= norm(v)
        vecs[:, r+1] .= v
    end
    return vals, vecs, m, (converged=converged, shift=shift)
end

# safe-shift wrapper: tries `shift0`, and on PosDefException reports back so the
# caller can re-assemble with a lower shift (H is destroyed by the attempt).
function ground_si_try(H::Matrix{Float64}, shift::Float64; kwargs...)
    try
        return ground_si(H, shift; kwargs...), :ok
    catch e
        e isa PosDefException || rethrow()
        return nothing, :not_spd
    end
end

end # module
