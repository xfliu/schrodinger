# run_checks.jl -- pre-sweep validation of the multi-centre sweep code path.
#
#   julia -t <threads> run_checks.jl
#
# Three checks, all on the H2+ system so recorded reference values exist:
#  A. PUBLISHED REPRODUCTION through the new code path at the published boxes and
#     their own npanel defaults, which the npanel rule round(4.8*Lx) reproduces
#     (48 at Lx=10, 96 at Lx=20):
#        Omega1 L=10 N=48 : mu1N_Ritz  vs [-0.5517637904077528,-0.5517637904035132]
#                           lambda1D   vs -0.550677376744471
#        Omega2 L=20 N=64 : lambda1D   vs -0.5509672617594417
#  B. KRYLOV vs LAPACK: shift-invert Lanczos ground eigenvalue against
#     eigen(Symmetric(A),1:1) on a moderate cell (D=15625), same matrix.
#  C. t_star SENSITIVITY: the same cell assembled at t_star = 1.0 and 1.5, to
#     confirm that raising t_star (the small-clearance guard) moves the
#     eigenvalue by less than the tail bound it removes.
using LinearAlgebra, Printf
include("sweep_core.jl")
using .SweepCore
const SC = SweepCore
const FC = SweepCore.FloatCore

BLAS.set_num_threads(haskey(ENV,"SWEEP_BLAS") ? parse(Int,ENV["SWEEP_BLAS"]) : 48)
@printf("threads=%d BLAS=%d julia=%s\n", Threads.nthreads(), BLAS.get_num_threads(), VERSION)
nodes, wts = FC.gl_rule(24)
cen2 = SC.h2plus_centres(2.0)
SHIFT = -0.5513170 - 0.35

println("\n=== A. published reproduction through the sweep code path ===")
for (L, N, what) in [(10.0, 48, "both"), (20.0, 64, "dirichlet")]
    LX = L; LY = 0.8L; LZ = 0.8L
    npan = SC.npanel_rule(LX)
    ts, tail = SC.choose_t_star(cen2, LX, LY, LZ; tol=1e-13)
    @printf("  L=%g N=%d npanel=%d (published default %s) t_star=%.2f tail=%.2e\n",
            LX, N, npan, L == 10 ? "48" : "96", ts, tail); flush(stdout)
    if what == "both"
        t1 = @elapsed (H,D) = SC.assemble_neumann_mc(N,:even,:even,:even; centres=cen2,
                LX=LX,LY=LY,LZ=LZ,npanel=npan,t_star=ts,nodes=nodes,wts=wts)
        t2 = @elapsed (v,_,it,inf1) = SC.ground_si(H, SHIFT; nev=2)
        H = nothing; GC.gc()
        @printf("    mu1N_Ritz = %.13f   (recorded interval [-0.5517637904077528,-0.5517637904035132])\n", v[1])
        @printf("    eee2      = %.10f   asm %.1fs eig %.1fs it=%d conv=%s\n", v[2], t1, t2, it, inf1.converged)
        flush(stdout)
    end
    t3 = @elapsed (HD,DD,nx) = SC.assemble_dirichlet_mc(N; centres=cen2,
            LX=LX,LY=LY,LZ=LZ,npanel=npan,t_star=ts,nodes=nodes,wts=wts)
    t4 = @elapsed (vd,_,itd,infd) = SC.ground_si(HD, SHIFT; nev=1)
    HD = nothing; GC.gc()
    rec = L == 10 ? -0.550677376744471 : -0.5509672617594417
    @printf("    lambda1D  = %.13f   recorded %.13f   diff %.2e   asm %.1fs eig %.1fs it=%d conv=%s\n",
            vd[1], rec, abs(vd[1]-rec), t3, t4, itd, infd.converged); flush(stdout)
end

println("\n=== B. Krylov shift-invert vs LAPACK dsyevr (same matrix, D=15625) ===")
let LX = 10.0, LY = 8.0, LZ = 8.0, N = 48
    npan = SC.npanel_rule(LX)
    (H,D) = SC.assemble_neumann_mc(N,:even,:even,:even; centres=cen2,
            LX=LX,LY=LY,LZ=LZ,npanel=npan,t_star=1.0,nodes=nodes,wts=wts)
    Hcopy = copy(H)
    tL = @elapsed esl = eigen(Symmetric(Hcopy), 1:2)
    Hcopy = nothing; GC.gc()
    tK = @elapsed (vk,_,itk,infk) = SC.ground_si(H, SHIFT; nev=2)
    H = nothing; GC.gc()
    @printf("  LAPACK  lambda1 = %.15f  lambda2 = %.13f   (%.1f s)\n", esl.values[1], esl.values[2], tL)
    @printf("  Krylov  lambda1 = %.15f  lambda2 = %.13f   (%.1f s, %d iters, conv=%s)\n",
            vk[1], vk[2], tK, itk, infk.converged)
    @printf("  |diff|  lambda1 = %.2e   lambda2 = %.2e   speedup = %.1fx\n",
            abs(vk[1]-esl.values[1]), abs(vk[2]-esl.values[2]), tL/tK); flush(stdout)
end

println("\n=== C. t_star sensitivity (guard is loss-free) ===")
let N = 40
    for (L, d, sysname) in [(10.0, 2.0, "H2plus"), (8.0, 4.0, "H3plus_d4")]
        cen = sysname == "H2plus" ? SC.h2plus_centres(2.0) : SC.h3plus_centres(d)
        LX = L; LY = 0.8L; LZ = 0.8L
        npan = SC.npanel_rule(LX)
        tsr, tailr = SC.choose_t_star(cen, LX, LY, LZ; tol=1e-13)
        vals = Float64[]
        for ts in (1.0, tsr)
            (H,D) = SC.assemble_neumann_mc(N,:even,:even,:even; centres=cen,
                    LX=LX,LY=LY,LZ=LZ,npanel=npan,t_star=ts,nodes=nodes,wts=wts)
            sh = (sysname == "H2plus" ? -0.5513170 : -0.76207999734) - 0.35
            (v,_,_,_) = SC.ground_si(H, sh; nev=1)
            H = nothing; GC.gc()
            push!(vals, v[1])
        end
        @printf("  %s L=%g N=%d: clearance=%.2f tail(t*=1)=%.2e -> rule t*=%.2f tail=%.2e\n",
                sysname, LX, N, LX - maximum(abs(c[1][1]) for c in cen),
                SC.cell_tail(cen, LX, LY, LZ, 1.0), tsr, tailr)
        @printf("    mu1(t*=1.00) = %.13f   mu1(t*=%.2f) = %.13f   shift = %.2e\n",
                vals[1], tsr, vals[2], abs(vals[2]-vals[1])); flush(stdout)
    end
end

println("CHECKS_DONE")
