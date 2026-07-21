
using Printf, LinearAlgebra
include("h2plus_bounds.jl")
using .H2plusBounds
import .H2plusBounds: make_matmat_sector, build_t_grid, BlockOp, DiagPrec
using IterativeSolvers

function build_aux(N::Int; nt::Int=48)
    t, wt = build_t_grid(nt)
    matmat, Kd, D = make_matmat_sector(N, t, wt, :even,:even,:even)
    return matmat, Kd, D
end
function eta_once(matmat, Kd, D, eps, s0, tol, maxiter)
    auxmat(X) = matmat(X) .+ (eps-1.0).*(Kd .* X) .+ s0 .* X
    A = BlockOp(auxmat, D)
    Pinv = DiagPrec(1.0 ./ (eps.*Kd .+ 5.0))
    r = lobpcg(A, false, 1; P=Pinv, tol=tol, maxiter=maxiter)
    return sort(r.λ)[1] - s0
end
function eta_of(matmat, Kd, D, eps)
    # retry across shifts/tols to dodge LOBPCG Gram breakdowns
    for (s0,tol) in [(30.0,1e-9),(45.0,1e-8),(20.0,1e-8),(60.0,1e-8),(35.0,5e-9)]
        try
            return eta_once(matmat, Kd, D, eps, s0, tol, 1000)
        catch e
            @printf("  (retry eps=%.2f s0=%.0f: %s)\n", eps, s0, typeof(e)); flush(stdout)
        end
    end
    return NaN
end
N = 128
matmat, Kd, D = build_aux(N)
@printf("# N=%d D=%d\n", N, D); flush(stdout)
for eps in [0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.80]
    eta = eta_of(matmat, Kd, D, eps)
    @printf("eps=%.2f  eta=%.6f  C_eps=%.6f\n", eps, eta, -eta); flush(stdout)
end
println("ALLDONE"); flush(stdout)
