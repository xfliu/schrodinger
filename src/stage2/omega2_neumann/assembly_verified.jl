# assembly_verified.jl — Track A / Stage-2, step A2
# Interval assembly of the shifted Neumann Galerkin sector Hamiltonian
#     Ĥ = K + P + σ I     (Matrix{Interval{Float64}})
# K   = kinetic diagonal ν_m (interval π²)
# P   = Coulomb = -(2/√π) Σ_t w_t (Bx_t ⊗ By_t ⊗ Bz_t)   [interval]
#       Bx_t = Bx_t^{(-2)} + Bx_t^{(+2)} (shifted sum);  By,Bz unshifted (s=0)
# σ I = optional coercivity shift (added by caller for the aux/PD operator)
#
# 1-D moment matrix from A1 frequency moments:
#   B_t(n,m) = 0.5 Nn Nm (G[|n-m|] + G[n+m]),  Nn=1/√(2L) (n=0) else 1/√L.
# Sector restriction: keep even (0,2,4,…) or odd (1,3,5,…) mode indices per axis.
module AssemblyVerified

using LinearAlgebra
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
include("moments_verified.jl")
using .MomentsVerified

const IV = Interval{Float64}
const LX, LY, LZ = 20.0, 16.0, 16.0
const AX = 2.0

# outer t-grid (Gauss–Legendre on (0,∞) via u/(1-u) map) — Stage-1 build_t_grid
function build_t_grid(nt::Int)
    k = 1:(nt-1); β = collect(k ./ sqrt.(4.0 .* k.^2 .- 1.0))
    J = SymTridiagonal(zeros(nt), β); vals, vecs = eigen(J)
    x = vals; w = 2.0 .* (vecs[1, :] .^ 2)
    u = 0.5 .* (x .+ 1); wu = 0.5 .* w
    return u ./ (1 .- u), wu ./ (1 .- u).^2
end

# interval 1-D moment matrix over frequency moments G[0..2N] at fixed t
# nodes/wts = certified small-t GL rule (from MomentsVerified.gl_reference)
function moment_matrix_1d(N::Int, L::Float64, s::Float64, t::Float64,
                          nodes::Vector{IV}, wts::Vector{IV}; npanel::Int=96, t_star::Float64=1.0)
    n1 = N+1; fmax = 2N
    G = Vector{IV}(undef, fmax+1)
    for f in 0:fmax
        κ = f*pi/(2L)
        G[f+1] = MomentsVerified.moment_verified(κ, t, L, s; t_star=t_star,
                                                  nodes=nodes, wts=wts, npanel=npanel)
    end
    Nn = [ (n==0) ? interval(1.0)/sqrt(interval(2.0)*interval(L)) :
                    interval(1.0)/sqrt(interval(L)) for n in 0:N ]
    B = Matrix{IV}(undef, n1, n1)
    half = interval(0.5)
    @inbounds for n in 0:N, m in 0:N
        B[n+1, m+1] = half * Nn[n+1] * Nn[m+1] * (G[abs(n-m)+1] + G[n+m+1])
    end
    return B
end

# select 1-based indices of even (0,2,…) or odd (1,3,…) modes
sector_idx(N::Int, p::Symbol) = p === :even ? collect(1:2:(N+1)) : collect(2:2:(N+1))

# Assemble Coulomb interval matrix P and kinetic diagonal Kdiag (interval),
# for the given parity sector. H = P + diag(Kdiag); aux = P + kin_scale*diag(Kdiag).
function assemble_PK(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                     nt::Int=48, npanel::Int=96, t_star::Float64=1.0,
                     nodes::Vector{IV}, wts::Vector{IV})
    t, wt = build_t_grid(nt)
    Ix, Iy, Iz = sector_idx(N,px), sector_idx(N,py), sector_idx(N,pz)
    nx, ny, nz = length(Ix), length(Iy), length(Iz)
    D = nx*ny*nz
    pref = interval(2.0)/sqrt(interval(pi))
    P = fill(interval(0.0), D, D)
    for k in 1:nt
        Bxm = moment_matrix_1d(N, LX, -AX, t[k], nodes, wts; npanel=npanel, t_star=t_star)
        Bxp = moment_matrix_1d(N, LX, +AX, t[k], nodes, wts; npanel=npanel, t_star=t_star)
        Bx = (Bxm .+ Bxp)[Ix, Ix]
        By = moment_matrix_1d(N, LY, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iy, Iy]
        Bz = moment_matrix_1d(N, LZ, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iz, Iz]
        wk = -pref * interval(wt[k])
        # column-major (x fastest) ⇒ global index = kron(Bz, kron(By, Bx))
        P .+= wk .* kron(Bz, kron(By, Bx))
    end
    νx = (interval(pi)/(interval(2.0)*interval(LX)))^2
    νy = (interval(pi)/(interval(2.0)*interval(LY)))^2
    νz = (interval(pi)/(interval(2.0)*interval(LZ)))^2
    mx = Ix .- 1; my = Iy .- 1; mz = Iz .- 1
    Kdiag = Vector{IV}(undef, D)
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        idx = (cz-1)*ny*nx + (cy-1)*nx + cx
        Kdiag[idx] = νx*interval(mx[cx]^2) + νy*interval(my[cy]^2) + νz*interval(mz[cz]^2)
    end
    return P, Kdiag, D
end

# Build interval Hamiltonian  H = kin_scale*diag(Kdiag) + P  (+ sigma*I).
function form_H(P::Matrix{IV}, Kdiag::Vector{IV}; kin_scale=interval(1.0), sigma=interval(0.0))
    D = length(Kdiag)
    H = copy(P)
    ks = kin_scale isa IV ? kin_scale : interval(kin_scale)
    sg = sigma isa IV ? sigma : interval(sigma)
    @inbounds for i in 1:D
        H[i,i] = H[i,i] + ks*Kdiag[i] + sg
    end
    return H
end

# Memory-lean assembler: identical intervals to assemble_PK, but fuses
# P .+= wk .* kron(...) into an in-place broadcast + forced GC each iteration,
# removing the extra D x D temporary (peak ~2 matrices instead of ~3).
function assemble_PK_lean(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                     nt::Int=48, npanel::Int=96, t_star::Float64=1.0,
                     nodes::Vector{IV}, wts::Vector{IV})
    t, wt = build_t_grid(nt)
    Ix, Iy, Iz = sector_idx(N,px), sector_idx(N,py), sector_idx(N,pz)
    nx, ny, nz = length(Ix), length(Iy), length(Iz)
    D = nx*ny*nz
    pref = interval(2.0)/sqrt(interval(pi))
    P = fill(interval(0.0), D, D)
    for k in 1:nt
        Bxm = moment_matrix_1d(N, LX, -AX, t[k], nodes, wts; npanel=npanel, t_star=t_star)
        Bxp = moment_matrix_1d(N, LX, +AX, t[k], nodes, wts; npanel=npanel, t_star=t_star)
        Bx = (Bxm .+ Bxp)[Ix, Ix]
        By = moment_matrix_1d(N, LY, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iy, Iy]
        Bz = moment_matrix_1d(N, LZ, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iz, Iz]
        wk = -pref * interval(wt[k])
        KB = kron(Bz, kron(By, Bx))        # one D x D allocation
        @inbounds @. P = P + wk * KB       # fused in-place: no extra temporary
        KB = nothing
        GC.gc()                             # reclaim before next D x D allocation
    end
    νx = (interval(pi)/(interval(2.0)*interval(LX)))^2
    νy = (interval(pi)/(interval(2.0)*interval(LY)))^2
    νz = (interval(pi)/(interval(2.0)*interval(LZ)))^2
    mx = Ix .- 1; my = Iy .- 1; mz = Iz .- 1
    Kdiag = Vector{IV}(undef, D)
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        idx = (cz-1)*ny*nx + (cy-1)*nx + cx
        Kdiag[idx] = νx*interval(mx[cx]^2) + νy*interval(my[cy]^2) + νz*interval(mz[cz]^2)
    end
    return P, Kdiag, D
end

end # module
