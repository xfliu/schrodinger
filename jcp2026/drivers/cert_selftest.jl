# cert_selftest.jl -- HARD GATE for the certified runs.
#
#   julia -t <threads> cert_selftest.jl
#
# Interval arithmetic is only safe to parallelise if it never depends on a global
# rounding mode.  IntervalArithmetic.jl 1.x is documented to use correctly-rounded
# operations without changing the FPU mode, but "documented" is not "verified on
# this host with this version", and the whole affordability of Phase 2 rests on it
# (D3 in cert_core.jl).  So this gate CHECKS it, bitwise, and exits non-zero if it
# fails.  The launcher refuses to start any production configuration unless this
# script exits 0.
#
# Three checks:
#   T1  the fused accumulation, run SERIALLY, is bit-identical to the published
#       expression  P .+= wkZ .* kron(Bz, kron(By, Bx))  on both endpoints.
#       (This is the algebraic-equivalence check: same association order.)
#   T2  the fused accumulation, run THREADED, is bit-identical to the serial run.
#       (This is the thread-safety check.)
#   T3  the threaded interval matvec is bit-identical to the serial one.
# The matrices are built from REAL certified moment matrices, not random data, so
# the exponent range and the interval widths are representative of production.
using LinearAlgebra, Printf
using IntervalArithmetic: Interval, interval, inf, sup, diam
push!(LOAD_PATH, get(ENV, "VEIGS_SRC", normpath(joinpath(@__DIR__, "..",
        "h2plus-rigorous-bounds", "lib", "Veigs.jl", "src"))))
using Veigs
include("cert_core.jl")
using .CertCore
const CC = CertCore
const IV = Interval{Float64}

bitsame(a::IV, b::IV) = (reinterpret(UInt64, inf(a)) == reinterpret(UInt64, inf(b))) &&
                        (reinterpret(UInt64, sup(a)) == reinterpret(UInt64, sup(b)))
function compare(A::Matrix{IV}, B::Matrix{IV}, name::String)
    bad = 0; first_bad = nothing
    @inbounds for j in axes(A,2), i in axes(A,1)
        if !bitsame(A[i,j], B[i,j])
            bad += 1
            first_bad === nothing && (first_bad = (i,j,A[i,j],B[i,j]))
        end
    end
    if bad == 0
        @printf("  %-46s PASS  (%d entries bit-identical)\n", name, length(A))
    else
        @printf("  %-46s FAIL  %d/%d entries differ\n", name, bad, length(A))
        i,j,a,b = first_bad
        @printf("      first at (%d,%d): [%.17g,%.17g] vs [%.17g,%.17g]\n",
                i, j, inf(a), sup(a), inf(b), sup(b))
    end
    return bad == 0
end

@printf("threads=%d julia=%s\n", Threads.nthreads(), VERSION)
@printf("IntervalArithmetic version: %s\n",
        try string(pkgversion(parentmodule(Interval))) catch; "unknown" end)
flush(stdout)

nodes, wts = CC.MomentsVerified.gl_reference(24)
N = 10                                    # D = 6^3 = 216 -> 46656 entries per matrix
cen = CC.h3plus_centres(4.0)              # 3 centres, exercises the grouping too
LX, LY, LZ = 12.0, 9.6, 9.6
npan = CC.npanel_rule(LX)
ts, _ = CC.choose_t_star(cen, LX, LY, LZ)
t, wt = CC.build_t_grid(48)
Ix = CC.sector_idx(N, :even)
Bx0 = CC.cos_moment_matrix_1d(N, LX, -4.0, t[7], nodes, wts; npanel=npan, t_star=ts)
Bx0 = Bx0 .+ CC.cos_moment_matrix_1d(N, LX, 0.0, t[7], nodes, wts; npanel=npan, t_star=ts)
Bx0 = Bx0 .+ CC.cos_moment_matrix_1d(N, LX,  4.0, t[7], nodes, wts; npanel=npan, t_star=ts)
Bx = Bx0[Ix,Ix]
By = CC.cos_moment_matrix_1d(N, LY, 0.0, t[7], nodes, wts; npanel=npan, t_star=ts)[Ix,Ix]
Bz = CC.cos_moment_matrix_1d(N, LZ, 0.0, t[7], nodes, wts; npanel=npan, t_star=ts)[Ix,Ix]
nx = length(Ix); D = nx^3
wkZ = -(interval(2.0)/sqrt(interval(pi))) * interval(wt[7])
@printf("test matrices: N=%d D=%d, moment interval widths %.2e .. %.2e\n",
        N, D, minimum(diam.(Bx)), maximum(diam.(Bx))); flush(stdout)

println("\n=== T1/T2: fused interval Kronecker accumulation ===")
Ppub = fill(interval(0.0), D, D)
Ppub .+= wkZ .* kron(Bz, kron(By, Bx))                     # the published expression
Pser = CC.kron_accum_iv!(fill(interval(0.0), D, D), wkZ, Bx, By, Bz; threaded=false)
Pthr = CC.kron_accum_iv!(fill(interval(0.0), D, D), wkZ, Bx, By, Bz; threaded=true)
ok1 = compare(Ppub, Pser, "T1 serial fused == published kron form")
ok2 = compare(Pser, Pthr, "T2 threaded fused == serial fused")

println("\n=== T3: threaded interval matvec ===")
x = [interval(sin(0.7*i) / sqrt(D)) for i in 1:D]
yser = CC.imatvec(Pser, x; threaded=false)
ythr = CC.imatvec(Pser, x; threaded=true)
ok3 = all(bitsame(yser[i], ythr[i]) for i in 1:D)
@printf("  %-46s %s  (%d entries)\n", "T3 threaded imatvec == serial imatvec",
        ok3 ? "PASS" : "FAIL", D)

println()
if ok1 && ok2 && ok3
    println("SELFTEST_PASS -- interval arithmetic is thread-exact on this host; threaded assembly authorised.")
    exit(0)
else
    println("SELFTEST_FAIL -- threaded interval arithmetic is NOT bit-identical.")
    println("Production runs must be re-launched with CERT_THREADED=0 (serial assembly,")
    println("~5.87e-6 * D^2 s per sector: ~2.1 h at D=35937 and ~7.8 h at D=68921).")
    exit(3)
end
