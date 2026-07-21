# pipeline2b.jl — precision-instrumented, memory-lean interval pipeline.
# Changes vs pipeline2.jl:
#   • sectors processed SEQUENTIALLY (eee assembled+consumed+freed, then oee) → ~half peak RAM,
#     lets N=64 complete where the simultaneous-sector version died.
#   • sharp-L evaluated with FULL interval propagation of μ1,μ2,C_ε (not point μ),
#     so sup(L)-inf(L) is an honest enclosure width for the lower bound.
#   • full-precision endpoints recorded for μ1^N,μ2^N,C_ε,η,L1,L2.
#   • inf(L_k) remains the rigorous lower bound (formula monotone: dL/dm>0, dL/dC<0;
#     interval inf takes the worst case over all propagated inputs).
module Pipeline2b
using LinearAlgebra, SparseArrays, Printf
include("assembly_verified.jl")
using .AssemblyVerified
const AV = AssemblyVerified
using .AssemblyVerified.MomentsVerified
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "..", "..", "lib", "Veigs.jl", "src")))
using Veigs
import Veigs: lehmann_behnke
using Arpack

const IV = Interval{Float64}
const LX = 20.0
const U1 = -0.545
nu_star(N::Int) = (interval(N+1)*interval(pi)/(interval(2.0)*interval(LX)))^2

# sharp-L with full interval inputs. m1 = ground Ritz interval, mk = target interval.
function sharp_L_iv(m1::IV, mk::IV, N::Int, eps::Float64, Ceps::IV)
    ε = interval(eps); one_ = interval(1.0); σ = Ceps
    ms1 = m1 + σ
    inf(ms1) > 0 || return nothing
    etaV = sqrt(ε/(one_-ε) + Ceps/ms1)
    g = etaV/sqrt(ms1)
    Ch2 = (one_ + g*g)/((one_-ε)*nu_star(N))
    msk = mk + σ
    return msk/(one_ + Ch2*msk) - σ
end

function lb_ground(Hi::Matrix{IV}, D::Int; nev::Int=3, rho::Float64, sigma::Float64)
    Hm = mid.(Hi)
    vals, vecs = eigs(Symmetric(Hm); nev=nev, which=:SR, maxiter=3000)
    p = sortperm(real(vals)); vals = real(vals[p]); vecs = real(vecs[:,p])
    V = reshape(vecs[:,1], D, 1)
    Bi = interval.(sparse(1.0I, D, D))
    t = @elapsed lg = lehmann_behnke(Hi, Bi, vals, V, rho, sigma, 1.0, 1, 1; do_shift=true)
    Hm = nothing
    return lg[1], t
end

function run_N(N::Int; nt::Int=48, npanel::Int=96, eps_grid=[0.40,0.30],
               nodes::Vector{IV}, wts::Vector{IV})
    T0 = time()
    # ---------- EEE block ----------
    t_as1 = @elapsed (Peee,Keee,Deee) = AV.assemble_PK_lean(N,:even,:even,:even; nt=nt,npanel=npanel,nodes=nodes,wts=wts)
    @printf("  [eee assemble] Deee=%d (%.1fs)\n",Deee,t_as1); flush(stdout)
    Heee = AV.form_H(Peee,Keee)
    μ1, t1 = lb_ground(Heee, Deee; rho=-0.60, sigma=-0.30); Heee=nothing; GC.gc()
    @printf("  [ritz eee] μ1^N=[%.14f,%.14f] w=%.3e (%.1fs)\n",inf(μ1),sup(μ1),diam(μ1),t1); flush(stdout)
    Ceps = Dict{Float64,IV}(); ηrec = Dict{Float64,Vector{Float64}}(); t_eta=0.0
    for e in eps_grid
        Haux = AV.form_H(Peee,Keee; kin_scale=interval(e))
        ηencl, te = lb_ground(Haux, Deee; rho=-1.60, sigma=-0.60); Haux=nothing; GC.gc()
        t_eta += te
        Ceps[e] = interval(-sup(ηencl), -inf(ηencl))
        ηrec[e] = [inf(ηencl), sup(ηencl)]
        @printf("  [eta eee] ε=%.2f η=[%.14f,%.14f] w=%.3e ⇒ C_ε=[%.14f,%.14f] (%.1fs)\n",
                e,inf(ηencl),sup(ηencl),diam(ηencl),inf(Ceps[e]),sup(Ceps[e]),te); flush(stdout)
    end
    Peee=nothing; Keee=nothing; GC.gc()
    # ---------- OEE block ----------
    t_as2 = @elapsed (Poee,Koee,Doee) = AV.assemble_PK_lean(N,:odd,:even,:even; nt=nt,npanel=npanel,nodes=nodes,wts=wts)
    @printf("  [oee assemble] Doee=%d (%.1fs)\n",Doee,t_as2); flush(stdout)
    Hoee = AV.form_H(Poee,Koee); Poee=nothing; Koee=nothing; GC.gc()
    μ2, t2 = lb_ground(Hoee, Doee; rho=-0.40, sigma=-0.25); Hoee=nothing; GC.gc()
    @printf("  [ritz oee] μ2^N=[%.14f,%.14f] w=%.3e (%.1fs)\n",inf(μ2),sup(μ2),diam(μ2),t2); flush(stdout)

    # ---------- sharp-L (full interval) ----------
    bL1=nothing;be1=NaN; bL2=nothing;be2=NaN
    for e in eps_grid
        L=sharp_L_iv(μ1,μ1,N,e,Ceps[e]); L===nothing && continue
        (bL1===nothing || inf(L)>inf(bL1)) && (bL1=L; be1=e)
    end
    for e in eps_grid
        L=sharp_L_iv(μ1,μ2,N,e,Ceps[e]); L===nothing && continue
        (bL2===nothing || inf(L)>inf(bL2)) && (bL2=L; be2=e)
    end
    @printf("  [sharpL] L1=[%.14f,%.14f] w=%.3e (ε=%.2f)\n",inf(bL1),sup(bL1),diam(bL1),be1); flush(stdout)
    @printf("  [sharpL] L2=[%.14f,%.14f] w=%.3e (ε=%.2f) sep=%s\n",inf(bL2),sup(bL2),diam(bL2),be2,inf(bL2)>U1); flush(stdout)
    br_lo=inf(bL1); br_hi=sup(μ1)
    @printf("  [bracket] λ1∈[%.14f,%.14f] w=%.3e (total %.1fs)\n",br_lo,br_hi,br_hi-br_lo,time()-T0); flush(stdout)

    return Dict{String,Any}(
        "N"=>N,"Deee"=>Deee,"Doee"=>Doee,
        "mu1N"=>[inf(μ1),sup(μ1)],"mu1N_w"=>diam(μ1),
        "mu2N"=>[inf(μ2),sup(μ2)],"mu2N_w"=>diam(μ2),
        "eta"=>Dict(string(e)=>ηrec[e] for e in eps_grid),
        "Ceps"=>Dict(string(e)=>[inf(Ceps[e]),sup(Ceps[e])] for e in eps_grid),
        "Ceps_w"=>Dict(string(e)=>diam(Ceps[e]) for e in eps_grid),
        "L1"=>[inf(bL1),sup(bL1)],"L1_w"=>diam(bL1),"eps1"=>be1,
        "L2"=>[inf(bL2),sup(bL2)],"L2_w"=>diam(bL2),"eps2"=>be2,
        "sep_L2_gt_U1"=>inf(bL2)>U1,
        "bracket_lo"=>br_lo,"bracket_hi"=>br_hi,"bracket_width"=>br_hi-br_lo,
        "t_assemble"=>t_as1+t_as2,"t_ritz"=>t1+t2,"t_eta"=>t_eta,"wall_total"=>time()-T0,
    )
end
end # module
