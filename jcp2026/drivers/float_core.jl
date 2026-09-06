# float_core.jl -- Float64 mirror of the Stage-2 interval pipeline.
#
# This is a DOUBLE-PRECISION counterpart of
#   src/stage2/<box>_neumann/{moments_verified.jl, assembly_verified.jl,
#                             dirichlet_verified.jl, lg_verified.jl}
# built for the discretisation-vs-truncation study.  It reproduces the MIDPOINTS
# of the interval quantities:
#   * the interval moment enclosure is  Ginf +- tail  (large t) or
#     compositeGL +- Bernstein  (small t); the float value is the centre term,
#     i.e. exactly what is computed here.
#   * the Neumann Galerkin mu1 / Dirichlet Ritz value / LG Moebius chain are
#     evaluated on the float matrices.
#
# Operator convention (identical to the interval code, verified against it):
#     H = -Delta - 1/|x-a1| - 1/|x-a2|,  a1,a2 = (-+AX,0,0),  AX = 2
#   -> cosine kinetic eigenvalue is (m*pi/(2L))^2 with coefficient ONE (-Delta,
#      not -(1/2)Delta), and the Coulomb prefactor is 2/sqrt(pi) with
#      coefficient ONE per centre.
#
# No external packages: LinearAlgebra + Printf only.
module FloatCore

using LinearAlgebra, Printf

# ---------------------------------------------------------------- GL rule ----
# Gauss-Legendre nodes/weights on [-1,1] by Golub-Welsch (Float64).
# Mirrors MomentsVerified.gl_reference(n) at double precision.
function gl_rule(n::Int)
    k = 1:(n-1)
    beta = collect(k ./ sqrt.(4.0 .* k .^ 2 .- 1.0))
    J = SymTridiagonal(zeros(n), beta)
    vals, vecs = eigen(J)
    return vals, 2.0 .* (vecs[1, :] .^ 2)
end

# outer t-grid: GL on (0,inf) via u/(1-u) map -- identical to
# AssemblyVerified.build_t_grid
function build_t_grid(nt::Int)
    k = 1:(nt-1)
    beta = collect(k ./ sqrt.(4.0 .* k .^ 2 .- 1.0))
    J = SymTridiagonal(zeros(nt), beta)
    vals, vecs = eigen(J)
    x = vals; w = 2.0 .* (vecs[1, :] .^ 2)
    u = 0.5 .* (x .+ 1); wu = 0.5 .* w
    return u ./ (1 .- u), wu ./ (1 .- u) .^ 2
end

# ------------------------------------------------------------- moments -------
# G(kappa,t) = int_{-L}^{L} cos(kappa (x+L)) exp(-t^2 (x-s)^2) dx
#
# large t (t >= t_star): infinite-domain closed form (the interval version adds
#                        a +-Gaussian-tail bound around exactly this value)
# small t             : composite Gauss-Legendre, npanel panels of half-length
#                       h = L/npanel, order n_p (the interval version adds a
#                       +-Bernstein remainder around exactly this value)
@inline function moment_larget(kappa::Float64, t::Float64, L::Float64, s::Float64)
    return cos(kappa * (s + L)) * (sqrt(pi) / t) * exp(-(kappa * kappa) / (4.0 * t * t))
end

function moment_smallt(kappa::Float64, t::Float64, L::Float64, s::Float64,
                       nodes::Vector{Float64}, wts::Vector{Float64}, npanel::Int)
    n_p = length(nodes)
    h = L / npanel
    acc = 0.0
    @inbounds for j in 0:(npanel - 1)
        cj = -L + (2j + 1) * h
        for i in 1:n_p
            x = cj + h * nodes[i]
            acc += (h * wts[i]) * cos(kappa * (x + L)) * exp(-t * t * (x - s) * (x - s))
        end
    end
    return acc
end

@inline function moment(kappa::Float64, t::Float64, L::Float64, s::Float64,
                        nodes::Vector{Float64}, wts::Vector{Float64};
                        npanel::Int, t_star::Float64=1.0)
    t >= t_star ? moment_larget(kappa, t, L, s) :
                  moment_smallt(kappa, t, L, s, nodes, wts, npanel)
end

# -------------------------------------------------- 1-D moment matrices ------
# Neumann cosine basis  phi_n(x) = N_n cos(n pi (x+L)/(2L)),  n = 0..N
#   N_0 = 1/sqrt(2L), N_n = 1/sqrt(L);   B[n,m] = 0.5 N_n N_m (G[|n-m|]+G[n+m])
function cos_moment_matrix_1d(N::Int, L::Float64, s::Float64, t::Float64,
                              nodes::Vector{Float64}, wts::Vector{Float64};
                              npanel::Int, t_star::Float64=1.0)
    n1 = N + 1; fmax = 2N
    G = Vector{Float64}(undef, fmax + 1)
    @inbounds for f in 0:fmax
        G[f + 1] = moment(f * pi / (2L), t, L, s, nodes, wts; npanel=npanel, t_star=t_star)
    end
    Nn = [(n == 0) ? 1.0 / sqrt(2.0 * L) : 1.0 / sqrt(L) for n in 0:N]
    B = Matrix{Float64}(undef, n1, n1)
    @inbounds for n in 0:N, m in 0:N
        B[n + 1, m + 1] = 0.5 * Nn[n + 1] * Nn[m + 1] * (G[abs(n - m) + 1] + G[n + m + 1])
    end
    return B
end

# Dirichlet sine basis  psi_n(x) = (1/sqrt(L)) sin(n pi (x+L)/(2L)),  n = 1..mmax
#   B[n,m] = 0.5 (1/L) (G[|n-m|] - G[n+m])      (product-to-sum: MINUS)
function sin_moment_matrix_1d(mmax::Int, L::Float64, s::Float64, t::Float64,
                              nodes::Vector{Float64}, wts::Vector{Float64};
                              npanel::Int, t_star::Float64=1.0)
    fmax = 2 * mmax
    G = Vector{Float64}(undef, fmax + 1)
    @inbounds for f in 0:fmax
        G[f + 1] = moment(f * pi / (2L), t, L, s, nodes, wts; npanel=npanel, t_star=t_star)
    end
    Ninv = 1.0 / L
    B = Matrix{Float64}(undef, mmax, mmax)
    @inbounds for n in 1:mmax, m in 1:mmax
        B[n, m] = 0.5 * Ninv * (G[abs(n - m) + 1] - G[n + m + 1])
    end
    return B
end

sector_idx(N::Int, p::Symbol) = p === :even ? collect(1:2:(N + 1)) : collect(2:2:(N + 1))
odd_modes(mmax::Int) = collect(1:2:mmax)

# ------------------------------------------ fused Kronecker accumulation -----
# P += w * kron(Bz, kron(By, Bx))  without allocating the D x D temporary.
# Column-major global index (x fastest): idx = ((cz-1)*ny + (cy-1))*nx + cx,
# matching AssemblyVerified's kron(Bz, kron(By, Bx)).
function kron_accum!(P::Matrix{Float64}, w::Float64,
                     Bx::Matrix{Float64}, By::Matrix{Float64}, Bz::Matrix{Float64})
    nx = size(Bx, 1); ny = size(By, 1); nz = size(Bz, 1)
    Threads.@threads for jz in 1:nz
        @inbounds for jy in 1:ny, jx in 1:nx
            jcol = ((jz - 1) * ny + (jy - 1)) * nx + jx
            for iz in 1:nz
                bz = w * Bz[iz, jz]
                bz == 0.0 && continue
                for iy in 1:ny
                    bzy = bz * By[iy, jy]
                    base = ((iz - 1) * ny + (iy - 1)) * nx
                    @simd for ix in 1:nx
                        P[base + ix, jcol] += bzy * Bx[ix, jx]
                    end
                end
            end
        end
    end
    return P
end

# ------------------------------------------------- Neumann eee assembly ------
# Returns the float sector Hamiltonian H = P + diag(K) (+ sigma I).
function assemble_neumann(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                          LX::Float64, LY::Float64, LZ::Float64, AX::Float64,
                          nt::Int=48, npanel::Int, t_star::Float64=1.0,
                          nodes::Vector{Float64}, wts::Vector{Float64},
                          kin_scale::Float64=1.0, sigma::Float64=0.0,
                          verbose::Bool=true)
    t, wt = build_t_grid(nt)
    Ix, Iy, Iz = sector_idx(N, px), sector_idx(N, py), sector_idx(N, pz)
    nx, ny, nz = length(Ix), length(Iy), length(Iz)
    D = nx * ny * nz
    pref = 2.0 / sqrt(pi)
    H = zeros(Float64, D, D)
    for k in 1:nt
        Bxm = cos_moment_matrix_1d(N, LX, -AX, t[k], nodes, wts; npanel=npanel, t_star=t_star)
        Bxp = cos_moment_matrix_1d(N, LX, +AX, t[k], nodes, wts; npanel=npanel, t_star=t_star)
        Bx = (Bxm .+ Bxp)[Ix, Ix]
        By = cos_moment_matrix_1d(N, LY, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iy, Iy]
        Bz = cos_moment_matrix_1d(N, LZ, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iz, Iz]
        kron_accum!(H, -pref * wt[k], Bx, By, Bz)
        if verbose && (k % 12 == 0)
            @printf("      [assemble] t-node %d/%d\n", k, nt); flush(stdout)
        end
    end
    nux = (pi / (2.0 * LX))^2; nuy = (pi / (2.0 * LY))^2; nuz = (pi / (2.0 * LZ))^2
    mx = Ix .- 1; my = Iy .- 1; mz = Iz .- 1
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        idx = ((cz - 1) * ny + (cy - 1)) * nx + cx
        H[idx, idx] += kin_scale * (nux * mx[cx]^2 + nuy * my[cy]^2 + nuz * mz[cz]^2) + sigma
    end
    return H, D
end

# ----------------------------------------------- Dirichlet eee assembly ------
# mmax = N+1 sine modes; even-in-x sector = ODD sine modes -> (N/2+1) per axis,
# matching the Neumann eee dimension.
function assemble_dirichlet(N::Int;
                            LX::Float64, LY::Float64, LZ::Float64, AX::Float64,
                            nt::Int=48, npanel::Int, t_star::Float64=1.0,
                            nodes::Vector{Float64}, wts::Vector{Float64},
                            verbose::Bool=true)
    mmax = N + 1
    t, wt = build_t_grid(nt)
    Ix = odd_modes(mmax); Iy = Ix; Iz = Ix
    nx = length(Ix); D = nx^3
    pref = 2.0 / sqrt(pi)
    H = zeros(Float64, D, D)
    for k in 1:nt
        Bxm = sin_moment_matrix_1d(mmax, LX, -AX, t[k], nodes, wts; npanel=npanel, t_star=t_star)
        Bxp = sin_moment_matrix_1d(mmax, LX, +AX, t[k], nodes, wts; npanel=npanel, t_star=t_star)
        Bx = (Bxm .+ Bxp)[Ix, Ix]
        By = sin_moment_matrix_1d(mmax, LY, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iy, Iy]
        Bz = sin_moment_matrix_1d(mmax, LZ, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star)[Iz, Iz]
        kron_accum!(H, -pref * wt[k], Bx, By, Bz)
        if verbose && (k % 12 == 0)
            @printf("      [assemble] t-node %d/%d\n", k, nt); flush(stdout)
        end
    end
    nux = (pi / (2.0 * LX))^2; nuy = (pi / (2.0 * LY))^2; nuz = (pi / (2.0 * LZ))^2
    @inbounds for cz in 1:nx, cy in 1:nx, cx in 1:nx
        idx = ((cz - 1) * nx + (cy - 1)) * nx + cx
        H[idx, idx] += nux * Ix[cx]^2 + nuy * Iy[cy]^2 + nuz * Iz[cz]^2
    end
    return H, D, nx
end

# --------------------------------------------------------- ground eigenpair --
# lowest eigenvalue + eigenvector of a dense symmetric matrix (LAPACK syevr
# restricted to index range 1:1 -- same eigenpair the reference obtains from a
# full eigen(Symmetric(.)), at a fraction of the back-transform cost).
function ground_pair(H::Matrix{Float64})
    es = eigen(Symmetric(H), 1:1)
    return es.values[1], vec(es.vectors[:, 1])
end

# ----------------------------------------------------------------- CG --------
# Conjugate gradients for the SPD dense system H w = b (mirrors the
# IterativeSolvers.cg! call in lg_verified.jl / run_lg_*.jl).
function cg_solve(H::Matrix{Float64}, b::Vector{Float64};
                  reltol::Float64=1e-12, maxiter::Int=8000)
    n = length(b)
    x = zeros(Float64, n)
    r = copy(b)
    p = copy(r)
    Ap = similar(r)
    rs = dot(r, r)
    nb = sqrt(rs)
    it = 0
    for i in 1:maxiter
        it = i
        mul!(Ap, Symmetric(H), p)
        alpha = rs / dot(p, Ap)
        axpy!(alpha, p, x)
        axpy!(-alpha, Ap, r)
        rs_new = dot(r, r)
        sqrt(rs_new) <= reltol * nb && break
        p .= r .+ (rs_new / rs) .* p
        rs = rs_new
    end
    return x, it
end

# --------------------------------------------------- LG Moebius (float) ------
# Single-test-vector Lehmann-Goerisch, A1 = <v,v> = 1:
#   A0 = <v,H'v>,  A2 = <v,H'^{-1}v>,  rho' = rho + c
#   A = A0 - rho',  B = A0 - 2 rho' + rho'^2 A2,  nu = A/B
#   L1_LG = rho' - rho'/(1-nu) - c
function lg_mobius(A0::Float64, A2::Float64, rho::Float64, c::Float64)
    rho_p = rho + c
    A = A0 - rho_p
    B = A0 - 2.0 * rho_p + rho_p * rho_p * A2
    nu = A / B
    lam_hat = rho_p - rho_p / (1.0 - nu)
    return lam_hat - c, A, B, nu
end

end # module
