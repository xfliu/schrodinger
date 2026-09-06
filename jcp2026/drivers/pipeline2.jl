# pipeline2.jl — Track A / Stage-2 final certified driver (A3–A5)
# Architecture validated at N=32:
#   • assemble interval P,K (kron)                          [A2]
#   • float ground vectors + 2nd-eee eigenvalue (Arpack)    [preconditioner only]
#   • lehmann_behnke DIRECT on H_eee, H_oee  -> verified μ1^N, μ2^N (rigorous)   [A3-Ritz]
#   • lehmann_behnke DIRECT on aux ε K + P   -> verified LOWER bound on η
#         => rigorous C_ε = -η   (σ = C_ε in sharp-L)                            [A3-η]
#   • sharp-L Möbius bound (interval arithmetic) -> rigorous L1, L2 (lower)       [A3-L]
#   • certified bracket  λ1 ∈ [ L1 , sup(μ1^N) ]                                  [A4/A5]
#
# Rigor directions:
#   μ_k^N (Ritz) : lehmann_behnke returns a verified ENCLOSURE of the k-th Galerkin
#       eigenvalue; sup(enclosure) is a rigorous UPPER bound on the true λ_k
#       (min-max), inf(enclosure) a rigorous lower bound on μ_k^N (Ritz value).
#   For sharp-L we need rigorous LOWER bounds on μ_k^N as the inputs (they play the
#       role of the projected eigenvalues); we use inf(enclosure).
#   η (aux ground energy): we need a rigorous LOWER bound on η so that
#       C_ε = -η is a rigorous UPPER bound on -η (safe, since dL/dC_ε < 0 ⇒ using a
#       larger C_ε only lowers L, keeping L a valid lower bound). lehmann_behnke's
#       inf(enclosure) is that rigorous lower bound on η.
module Pipeline2

using LinearAlgebra, SparseArrays, Printf
include("assembly_verified.jl")
using .AssemblyVerified
const AV = AssemblyVerified
using .AssemblyVerified.MomentsVerified
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
push!(LOAD_PATH, expanduser("~/Workspace/Julia_LIB/Lib/Veigs.jl/src"))
using Veigs
import Veigs: lehmann_behnke
using Arpack

const IV = Interval{Float64}
const LX = 10.0
const U1 = -0.545
const LAM1 = -0.551317
const LAM2 = -0.333768

nu_star(N::Int) = (interval(N+1)*interval(pi)/(interval(2.0)*interval(LX)))^2

# sharp-L Möbius (interval). μ1_lo, μk_lo = rigorous lower bounds on Ritz eigenvalues.
# Ceps = rigorous upper bound (interval) on C_ε (= -η). ε exact rational-ish float.
function sharp_L(μ1_lo::Float64, μk_lo::Float64, N::Int, eps::Float64, Ceps::IV)
    ε = interval(eps); one_ = interval(1.0)
    σ = Ceps
    ms1 = interval(μ1_lo) + σ
    inf(ms1) > 0 || return nothing
    etaV = sqrt(ε/(one_-ε) + Ceps/ms1)
    g = etaV/sqrt(ms1)
    Ch2 = (one_ + g*g)/((one_-ε)*nu_star(N))
    msk = interval(μk_lo) + σ
    Lk = msk/(one_ + Ch2*msk) - σ
    return Lk    # inf(Lk) is a rigorous lower bound on λ_k
end

# verified enclosure of the ground (r=s=1) eigenvalue of interval matrix Hi (identity mass)
# using lehmann_behnke seeded by float spectrum. Returns (enclosure::IV, wall).
function lb_ground(Hi::Matrix{IV}, Deee::Int; nev::Int=3, rho::Float64, sigma::Float64)
    Hm = mid.(Hi)
    vals, vecs = eigs(Symmetric(Hm); nev=nev, which=:SR, maxiter=3000)
    p = sortperm(real(vals)); vals = real(vals[p]); vecs = real(vecs[:,p])
    V = reshape(vecs[:,1], Deee, 1)
    Bi = interval.(sparse(1.0I, Deee, Deee))
    t = @elapsed lg = lehmann_behnke(Hi, Bi, vals, V, rho, sigma, 1.0, 1, 1; do_shift=true)
    return lg[1], vals, t
end

function run_N(N::Int; nt::Int=48, npanel::Int=48,
               eps_grid::Vector{Float64}=[0.40,0.30],
               nodes::Vector{IV}, wts::Vector{IV}, verbose::Bool=true)
    T0 = time()
    t_as = @elapsed begin
        Peee,Keee,Deee = AV.assemble_PK(N,:even,:even,:even; nt=nt, npanel=npanel, nodes=nodes, wts=wts)
        Poee,Koee,Doee = AV.assemble_PK(N,:odd, :even,:even; nt=nt, npanel=npanel, nodes=nodes, wts=wts)
    end
    verbose && @printf("  [assemble] Deee=%d Doee=%d (%.1fs)\n",Deee,Doee,t_as); flush(stdout)

    Heee = AV.form_H(Peee,Keee); Hoee = AV.form_H(Poee,Koee)

    # A3-Ritz: verified μ1^N (eee), μ2^N (oee). Separators are generous rough barriers
    # (below λ1 for ρ; between λ1 and λ2^sector for σ). Physical λ1≈-0.5513, gaps large.
    μ1, v1, t1 = lb_ground(Heee, Deee; rho=-0.60, sigma=-0.30)
    μ2, v2, t2 = lb_ground(Hoee, Doee; rho=-0.40, sigma=-0.25)
    μ1_lo, μ1_hi = inf(μ1), sup(μ1)
    μ2_lo, μ2_hi = inf(μ2), sup(μ2)
    verbose && @printf("  [ritz] μ1^N=[%.10f,%.10f](%.1fs) μ2^N=[%.10f,%.10f](%.1fs)\n",
                       μ1_lo,μ1_hi,t1,μ2_lo,μ2_hi,t2); flush(stdout)

    # A3-η: verified lower bound on aux ground energy η for each ε; C_ε = -η (interval).
    # aux ground η ≈ -0.9..-1.1; 2nd aux eig ≈ -0.38..-0.51 ⇒ σ=-0.60 separates, ρ=-1.5 barrier.
    Ceps = Dict{Float64,IV}(); t_eta = 0.0
    for e in eps_grid
        Haux = AV.form_H(Peee,Keee; kin_scale=interval(e))
        ηencl, _, te = lb_ground(Haux, Deee; rho=-1.60, sigma=-0.60)
        t_eta += te
        # η ≥ inf(ηencl)  ⇒  -η ≤ -inf(ηencl); C_ε upper bound = -inf(ηencl).
        Ceps[e] = interval(-sup(ηencl), -inf(ηencl))   # rigorous interval for C_ε = -η
        verbose && @printf("  [eta] ε=%.2f η∈[%.8f,%.8f] ⇒ C_ε∈[%.6f,%.6f] (%.1fs)\n",
                           e, inf(ηencl), sup(ηencl), inf(Ceps[e]), sup(Ceps[e]), te); flush(stdout)
    end

    # A3-L: sharp-L lower bounds, optimize ε over grid (max inf).
    bestL1=nothing; be1=NaN
    for e in eps_grid
        L=sharp_L(μ1_lo, μ1_lo, N, e, Ceps[e]); L===nothing && continue
        if bestL1===nothing || inf(L)>inf(bestL1); bestL1=L; be1=e; end
    end
    bestL2=nothing; be2=NaN
    for e in eps_grid
        L=sharp_L(μ1_lo, μ2_lo, N, e, Ceps[e]); L===nothing && continue
        if bestL2===nothing || inf(L)>inf(bestL2); bestL2=L; be2=e; end
    end
    L1_lo = inf(bestL1); L2_lo = inf(bestL2)
    sep = L2_lo > U1
    verbose && @printf("  [sharpL] L1≥%.6f (ε=%.2f) L2≥%.6f (ε=%.2f) sep(L2>U1)=%s\n",
                       L1_lo,be1,L2_lo,be2,sep); flush(stdout)

    # Certified bracket for λ1:  [ L1_lo , μ1_hi ]  (μ1_hi rigorous upper bound on λ1).
    br_lo = L1_lo; br_hi = μ1_hi; br_w = br_hi - br_lo
    verbose && @printf("  [bracket] λ1 ∈ [%.8f, %.8f]  width=%.3e  (total %.1fs)\n",
                       br_lo,br_hi,br_w,time()-T0); flush(stdout)

    return Dict{String,Any}(
        "N"=>N,"Deee"=>Deee,"Doee"=>Doee,
        "mu1N_lo"=>μ1_lo,"mu1N_hi"=>μ1_hi,"mu2N_lo"=>μ2_lo,"mu2N_hi"=>μ2_hi,
        "Ceps"=>Dict(string(e)=>[inf(Ceps[e]),sup(Ceps[e])] for e in eps_grid),
        "L1_lo"=>L1_lo,"eps1"=>be1,"L2_lo"=>L2_lo,"eps2"=>be2,"sep_L2_gt_U1"=>sep,
        "bracket_lo"=>br_lo,"bracket_hi"=>br_hi,"bracket_width"=>br_w,
        "t_assemble"=>t_as,"t_ritz"=>t1+t2,"t_eta"=>t_eta,"wall_total"=>time()-T0,
    )
end

end # module
