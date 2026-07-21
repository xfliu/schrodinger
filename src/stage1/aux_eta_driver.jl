
using Printf, LinearAlgebra
include("h2plus_bounds.jl")
using .H2plusBounds
import .H2plusBounds: make_matmat_sector, build_t_grid, BlockOp, DiagPrec
using IterativeSolvers

# Sharp form constant via eta:  eta = inf [ eps*||grad u||^2 - (V^- u,u) ] / ||u||^2
# For H2+, V = -sum Z_i/|x-a_i| so V^- = -V (all attractive), (V^- u,u) = -(V u,u).
# Our matmat implements  ||grad||^2 (Kdiag) + (V u,u)   [Coulomb sign already negative].
# So  eps*||grad||^2 - (V^- u,u) = eps*Kdiag + (V u,u) = matmat(X) + (eps-1)*Kdiag.*X.
# eta = ground eigenvalue of that operator (Neumann, eee sector -> constant mode lives here).

function aux_eta(N::Int, eps::Float64; nt::Int=48, tol=1e-9, maxiter=600)
    t, wt = build_t_grid(nt)
    matmat, Kd, D = make_matmat_sector(N, t, wt, :even,:even,:even)
    s0 = 60.0   # PD shift for lobpcg: operator + s0 I  (eta expected ~ -10..-20)
    auxmat(X) = matmat(X) .+ (eps-1.0).*(Kd .* X) .+ s0 .* X
    A = BlockOp(auxmat, D)
    Pinv = DiagPrec(1.0 ./ (eps.*Kd .+ s0))
    r = lobpcg(A, false, 1; P=Pinv, tol=tol, maxiter=maxiter)
    return sort(r.λ)[1] - s0, D
end

# validate at a modest N and the eps used in the analytic optimum
for (N,eps) in [(64,0.408),(96,0.408),(128,0.408)]
    eta, D = aux_eta(N, eps)
    @printf("N=%d eps=%.3f  D=%d  eta=%.6f  C_eps=-eta=%.6f\n", N, eps, D, eta, -eta)
end
