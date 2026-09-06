# run_cert.jl -- ONE driver per configuration, ONE output file containing the
# final certified interval directly.
#
#   julia -t <threads> run_cert.jl <CONFIG> [--serial-check]
#
# CONFIG:
#   C1        H2+      L=12  N=48  N'=64
#   C2        H2+      L=14  N=64  N'=80
#   C3        H3^2+    d=4   L=12  N=64  N'=80
#   SELFTEST  H2+      L=10  N=32  N'=48   (the published Omega1 configuration,
#             run with the PUBLISHED lehmann_behnke rho/sigma so that Stage-A must
#             reproduce results/certified_precision.json and the LG must reproduce
#             the certified_lg_results.csv row for primary_N=32, aux=48)
#
# Every stage checkpoints its own JSON as it completes, so a late failure never
# loses an earlier stage.  The final file cert_<CONFIG>.json carries the complete
# chain AND the final enclosure in one place -- no cross-file combination and no
# re-evaluated Moebius step is needed to read off the headline interval.
using LinearAlgebra, Printf, SparseArrays, JSON
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
using Arpack, IterativeSolvers
# Veigs is a path package in the bundle; VEIGS_SRC is set by the launcher.
push!(LOAD_PATH, get(ENV, "VEIGS_SRC", normpath(joinpath(@__DIR__, "..",
        "h2plus-rigorous-bounds", "lib", "Veigs.jl", "src"))))
using Veigs
include("sigma_loc_iv.jl")
include("cert_core.jl")
include("lg_goerisch_correct2.jl")
using .SigmaLoc
using .CertCore
using .GoerischCorrect2
const CC = CertCore
const G2 = GoerischCorrect2
const IV = Interval{Float64}

# ===================== CORRECTED LEHMANN-GOERISCH REALIZATION =====================
# The shipped realization takes X = V_{N'}, T = zero-padding, b_G = a_c on V_{N'}.  (A3)
# then holds only for u,v in V_{N'} and (A4) is enforced only against v in V_{N'}, so the
# theorem applies with D = V_{N'} and bounds mu_{k,N'} + c from below -- not the operator's
# eigenvalues.  This driver replaces the Goerisch space by
#     X = (L^2)^3 x L^2 x H^1(Omega),  Tu = {sqrt(1-eps') grad u, sqrt(c-C_eps') u, u},
#     b_G = (p,p') + (s,s') + [eps'(grad q,grad q') + ((V+C_eps')q,q')],
# which keeps D = H^1(Omega) and is >= 0 on ALL of X.  See lg_goerisch_correct2.jl.
#
# ENV CONTROLS
#   CERT_G_CEPSFILE  run_ceps.jl output supplying certified C_eps' per eps'  (REQUIRED)
#   CERT_G_EPS       comma list of eps' to evaluate (default: every eps' in that file)
#   CERT_G_DISC      "0" forces the DISCRETIZATION term ||(I-Pi^0)(V wtilde)||^2 to zero.
#                    This is validation gate G1 only: with it the corrected code path must
#                    reproduce the recorded certificates bit-identically, which proves the
#                    new code changes nothing except by adding the term it is meant to add.
#   CERT_G_FROM      prefix of an EXISTING certificate's stage files (e.g. /path/cert_SL1).
#                    Stage A and the Dirichlet stage are unchanged by this repair -- neither
#                    uses Lehmann-Goerisch -- so they are read from that configuration's OWN
#                    recorded stage files instead of recomputed, with the geometry, N and
#                    Ceps grid asserted to match.  Unset = recompute everything.
#   CERT_G_OFF       "1" runs the shipped uncorrected A_2 (for A/B comparison).
const GCEPSFILE = get(ENV, "CERT_G_CEPSFILE", "")
const GEPSLIST  = haskey(ENV,"CERT_G_EPS") ?
    [parse(Float64,x) for x in split(ENV["CERT_G_EPS"], ",")] : Float64[]
const GDISC     = get(ENV, "CERT_G_DISC", "1") != "0"
const GFROM     = get(ENV, "CERT_G_FROM", "")
const GOFF      = get(ENV, "CERT_G_OFF", "0") == "1"

# certified C_eps' per eps', read as a POINT interval at the certified upper bound so no
# value below C_eps^opt can be consumed.  The file's geometry is checked against the
# configuration in gate_G4! below.
const GCEPS = Dict{Float64,IV}()
const GCEPSREC = Dict{String,Any}()
function load_gceps!()
    GCEPSFILE == "" && return
    J = JSON.parsefile(GCEPSFILE)
    for (ek, ev) in J["results"]
        # JSON.parsefile here yields JSON.Object, NOT Dict, for nested objects: an
        # `ev isa Dict` guard silently skips every entry (this is exactly what defeated
        # the earlier goer/run_goerisch.jl run, which reported "no certified Ceps").
        (ev isa AbstractDict && haskey(ev, "Ceps_cert")) || continue
        GCEPS[parse(Float64, String(ek))] = interval(Float64(ev["Ceps_cert"]))
        GCEPSREC[String(ek)] = ev
    end
    GCEPSREC["_file"] = GCEPSFILE
    GCEPSREC["_box"] = J["box"]; GCEPSREC["_centres"] = J["centres"]
    GCEPSREC["_N_aux"] = J["N_aux"]; GCEPSREC["_system"] = J["system"]
    GCEPSREC["_d_paper"] = J["d_paper"]
    return
end

# G4 SELF-CONSISTENCY, machine-checkable: every quantity entering A_2^G must come from this
# configuration's own (box, centres, N', sector, sigma).  Anything read from a file is
# checked field by field against the live configuration and the booleans are recorded.
const G4 = Dict{String,Any}()
function gate_G4!(key::String, checks::Vector{Pair{String,Bool}}, values::Dict{String,Any})
    d = Dict{String,Any}("checks"=>Dict(k=>v for (k,v) in checks),
                         "values"=>values, "all_pass"=>all(v for (_,v) in checks))
    G4[key] = d
    for (k,v) in checks
        v || @printf("  [G4 FAIL] %s : %s\n", key, k)
    end
    @assert d["all_pass"] "G4 self-consistency failed for $key: $(d)"
    return d
end

eps_grid_for_goerisch() = isempty(GEPSLIST) ? sort(collect(keys(GCEPS)); rev=true) : GEPSLIST

# stage_A's own resolution rule, exposed so a from-file read can assert like with like
NSTAGEA(c) = haskey(ENV,"CERT_NSTAGEA") ? parse(Int, ENV["CERT_NSTAGEA"]) : c.N_stageA

"""
    goerisch_block(...) -> (rec, cands)

Everything the corrected realization needs, for ONE stage.  `cands` is the list of
(eps', Ceps', A2^G) over the admissible eps'; every entry is a VALID certificate, so the
caller may take the sharpest exactly as the shipped driver already takes the sharpest
separator and the sharpest sharp-L eps.
"""
function goerisch_block(label::String, wI::Vector{IV}, Hw::Vector{IV}, r::Vector{IV},
                        Kdiag::Vector{IV}, wv::IV, wHw::IV, r2sq::IV,
                        cshift::Float64, Nprime::Int, LX::Float64,
                        nx::Int, ny::Int, nz::Int, momx::Function, momyz::Function,
                        centres::Vector{CC.Centre}, tgrid::Vector{Float64},
                        twts::Vector{Float64}, skip_const::Bool,
                        prov::Dict{String,Any})
    D = length(wI)
    # Pi^0_{N'}(V wtilde) coefficients = H'w - (Kdiag + c)w, free given Hw.  Valid enclosure:
    # form_H_inplace! only ever added ks*Kdiag[i] + sigma to the DIAGONAL.
    Pw = Vector{IV}(undef, D)
    @inbounds for i in 1:D; Pw[i] = Hw[i] - (Kdiag[i] + interval(cshift))*wI[i]; end
    projsq = CC.idot(Pw, Pw)
    r0 = skip_const ? r[1] : interval(0.0)
    Pw = nothing
    @printf("  [%s G] ||r||^2 <= %.6e  r_0 = %+.3e  ||Pi^0(V w)||^2 = [%.14f,%.14f]\n",
            label, sup(r2sq), mid(r0), inf(projsq), sup(projsq)); flush(stdout)

    tv = 0.0; vsq = interval(0.0); nterms = 0; npairs = 0
    if GDISC
        tv = @elapsed ((vsq, nterms, npairs) =
            G2.vsq_quadform(wI, nx, ny, nz, momx, momyz, centres, tgrid, twts))
        tailsq = vsq - projsq
        @printf("  [%s G] ||V w||^2 = [%.14f,%.14f] (%d terms / %d node pairs, %.1f s)\n",
                label, inf(vsq), sup(vsq), nterms, npairs, tv)
    else
        tailsq = interval(0.0)     # GATE G1: discretization term forced to zero
        @printf("  [%s G] CERT_G_DISC=0 -- discretization term FORCED TO ZERO (gate G1)\n", label)
    end
    flush(stdout)

    nus  = G2.nu_star(Nprime, LX)
    nu1n = G2.nu1_neumann(LX)
    gphi, inpart, outpart, ts = G2.grad_phi_sq(r, Kdiag, tailsq, nus; skip_const=skip_const)
    dom = sup(outpart) > sup(inpart) ? "discretization" : "algebraic"
    @printf("  [%s G] |grad phi|^2 <= %.6e  = algebraic %.6e + discretization %.6e (%s dominates)\n",
            label, sup(gphi), sup(inpart), sup(outpart), dom)
    @printf("        nu_*(N'=%d) = %.8f (driver convention ((N'+1)pi/(2LX))^2) ; nu_1 = %.8f\n",
            Nprime, mid(nus), mid(nu1n)); flush(stdout)

    cands = Any[]
    perreps = Dict{String,Any}()
    for e in eps_grid_for_goerisch()
        ce = get(GCEPS, e, nothing)
        adm = ce !== nothing && sup(ce) < cshift
        if !adm
            perreps[string(e)] = Dict{String,Any}("admissible"=>false,
                "Ceps_cert"=>ce === nothing ? nothing : sup(ce),
                "c"=>cshift,
                "reason"=>ce === nothing ? "no certified Ceps' at this eps'" :
                          "c = $(cshift) does not exceed Ceps' = $(sup(ce))")
            @printf("  [%s G] eps'=%.4g : NOT ADMISSIBLE (Ceps' = %s, c = %.4f)\n", label, e,
                    ce === nothing ? "missing" : string(sup(ce)), cshift); flush(stdout)
            continue
        end
        corr, defect, sigma2 = G2.a2_correction(gphi, e, cshift, ce, r0; has_const=skip_const)
        A2G = interval(2.0)*wv - wHw + corr
        push!(cands, (e, ce, A2G, corr, defect, sigma2))
        perreps[string(e)] = Dict{String,Any}("admissible"=>true, "Ceps_cert"=>sup(ce),
            "c"=>cshift, "c_minus_Ceps"=>cshift - sup(ce),
            "A2_goerisch"=>[inf(A2G), sup(A2G)], "correction_sup"=>sup(corr),
            "defect_sup"=>sup(defect), "sigma2_const_mode_sup"=>sup(sigma2))
        @printf("  [%s G] eps'=%.4g Ceps'=%.8f (c-Ceps'=%+.6f) : A2^G=[%.12f,%.12f] corr=%.6e\n",
                label, e, sup(ce), cshift - sup(ce), inf(A2G), sup(A2G), sup(corr)); flush(stdout)
    end

    rec = Dict{String,Any}(
      "realization"=>"corrected: X = (L^2)^3 x L^2 x H^1(Omega), D = H^1(Omega)",
      "discretization_term_included"=>GDISC,
      "Vw_norm_sq"=>GDISC ? [inf(vsq),sup(vsq)] : nothing,
      "proj_norm_sq"=>[inf(projsq),sup(projsq)],
      "tail_sq"=>[inf(ts),sup(ts)],
      "nu_star"=>mid(nus), "nu_star_convention"=>"((N'+1)*pi/(2*LX))^2, the driver's own conservative value",
      "nu_1_neumann"=>mid(nu1n),
      "nu_1_note"=>("in-space part uses the EXACT modal sum sum_{k!=0} |r_k|^2/nu_k, which is " *
                    "sharper than and implies the brief's ||r||^2/nu_1; nu_1 is reported for reference"),
      "grad_phi_sq_total_sup"=>sup(gphi),
      "grad_phi_sq_algebraic_sup"=>sup(inpart),
      "grad_phi_sq_discretization_sup"=>sup(outpart),
      "dominant_term"=>dom,
      "r0_const_mode"=>mid(r0), "r2_sup"=>sup(r2sq),
      "vsq_terms"=>nterms, "vsq_node_pairs"=>npairs, "t_vsq"=>tv,
      # parameters as SEEN INSIDE this block, recorded by the moment closures themselves at the
      # moment they were called.  gate_G4! compares these against the stage's own variables, so
      # the provenance check is a runtime comparison of two independently-sourced values rather
      # than an assertion that the call site passed what it passed.
      "Nprime_used"=>Nprime, "nx_used"=>nx, "ny_used"=>ny, "nz_used"=>nz,
      "npanel_used"=>get(prov,"npanel",nothing), "t_star_used"=>get(prov,"t_star",nothing),
      "LX_used"=>get(prov,"LX",nothing), "LYZ_used"=>get(prov,"LYZ",nothing),
      "nt_used"=>length(tgrid), "momx_calls"=>get(prov,"nx_calls",0),
      "momyz_calls"=>get(prov,"nyz_calls",0),
      "per_eps_prime"=>perreps)
    return rec, cands
end

# Goerisch shift c: H' = H + cI must be CERTIFIED positive definite, which requires
# c + inf(L1) > 0 where inf(L1) is the Stage-A rigorous lower bound on lambda_1.
# c = 1.0 is an H2+-TUNED value: for H2+ lambda_1 ~ -0.551 and inf(L1) ~ -0.61, so c=1
# leaves a margin of ~0.39. For linear H3^2+ at d_paper=2.8, lambda_1 ~ -0.9365 and the
# Stage-A bound at the separator resolution is inf(L1) = -1.0515, so c=1 gives
# c + inf(L1) = -0.0515 and the SPD certification FAILS. c is a free parameter of the
# Goerisch construction -- any c making H+cI certified positive definite is admissible --
# so it is set adaptively from the certified Stage-A bound rather than hardcoded.
# This is a prerequisite of the method, NOT a tuning knob for the resulting inequality.
const CSHIFT_REF = Ref(1.0)
CSHIFT() = CSHIFT_REF[]
function set_shift_from_L1!(L1::IV; margin::Float64=0.25)
    need = -inf(L1) + margin
    CSHIFT_REF[] = max(1.0, ceil(need*4)/4)      # quarter-integer, never below the H2+ value
    return CSHIFT_REF[]
end
# EPS_GRID may be overridden with CERT_EPS (comma separated). The driver takes the BEST
# (largest) sharp-L over the grid, so adding a value can only help; it is not a tuning knob
# on the result, and every grid member is reported with its own Ceps and sharp-L output.
const EPS_GRID = haskey(ENV,"CERT_EPS") ?
    [parse(Float64,x) for x in split(ENV["CERT_EPS"], ",")] : [0.40, 0.30]

# sigma_loc parameters per eps.  delta is given as an exact Rational so the interval encloses
# the intended rational value rather than the nearest Float64.  r1 = rho/2, r2 = rho, with rho
# computed from the geometry (box_rho).  On Omega_1 (rho=8) this is r1=4, r2=8, reproducing the
# author-specified parameters exactly; both (eps,delta) pairs happen to be the stationary point
# of sigma_loc in delta, so the constant is insensitive to delta at first order.
# The 0.50 entry is a FALLBACK only: it is consulted solely on the analytic sigma_loc path,
# which CERT_CEPS_FILE bypasses. Its value is the closed form delta* = 2 alpha/(Ztot (r2-r1))
# evaluated on Omega_1 (Ztot=2, r2-r1=4), so on any other geometry the analytic path would be
# valid but not delta-optimal -- acceptable for a fallback that the certified runs never take.
const SIGMA_LOC_DELTA = Dict(0.40 => 1//10, 0.30 => 3//40, 0.50 => 1//8)

# ---------------------------------------------------------------- configs ----
struct Config
    label::String; system::String; d::Float64
    L::Float64; N::Int; Nprime::Int
    N_sep::Int; Nprime_sep::Int      # resolution for the certified mu2 separator stage
    N_stageA::Int                    # resolution for Stage A.  DECOUPLED from N_sep: the new
                                     # sigma_loc shift needs N_stageA=80 for L_2 > U_1, but
                                     # raising the mu2 stage to match would force its oee
                                     # auxiliary to N'>=96 (dim 115248, 212 GB) -- impossible.
    lambda_ref::Float64; lambda_ref_err::Float64; ref_kind::String
    published_rs::Bool          # use the published lehmann_behnke rho/sigma
    float_width::Float64        # W at this cell from optimal_box_surface.csv
    float_lower::Float64        # LG_lower  at this cell (float sweep)
    float_upper::Float64        # lambda1D_Ritz at this cell (float sweep)
    ref_var_upper::Float64      # independently verified variational UPPER bound on lambda_1, or NaN
    purpose::String
end
function get_config(tag::String)
    tag == "C1" && return Config("C1","H2plus",2.0,12.0,48,64,32,48,32,
        -0.5513170, 0.0, "exact_prolate_spheroidal_half_of_-1.1026342", false,
        0.0004715610172789, -0.5515394845756576, -0.55106792355837864,
        NaN,
        "cost result: the published enclosure width reached at N=48/N'=64 instead of N=64/N'=80")
    tag == "C2" && return Config("C2","H2plus",2.0,14.0,64,80,32,48,32,
        -0.5513170, 0.0, "exact_prolate_spheroidal_half_of_-1.1026342", false,
        0.0002236010141072, -0.55139271073453955, -0.5511691097204322,
        NaN,
        "sharpness result: 2.13x tighter than the published 4.763392e-4 at identical N, N' and cost")
    # (3a) de-risking run: the three-centre INTERVAL assembly has not been run in
    # production before, only as a moment containment test and a zero-charge
    # regression.  Cheap, and it fails fast if anything is wrong.
    tag == "C3a" && return Config("C3a","H3plus",4.0,12.0,32,48,32,48,32,
        -0.76207999734, 1.5e-7, "floating_gaussian_variational_upper_bound", false,
        0.0013106719624234001, -0.76265656780735147, -0.76134589584492807,
        NaN,
        "de-risking: first production run of the three-centre interval assembly")
    # (3c) certificate ROBUSTNESS ACROSS N on the three-centre operator: a second,
    # coarser N so the B>0 / nu<1 conditions are shown to hold across resolutions
    # rather than at one cell.  This cell was never run at LG tier in the float
    # sweep, so it has NO float LG lower bound and NO float width -- those stay NaN
    # and the driver emits null rather than a fabricated comparison.  If the
    # certificate fails here, that is the reportable finding; nothing is tuned.
    tag == "C3c" && return Config("C3c","H3plus",4.0,12.0,24,40,24,40,24,
        -0.76207999734, 1.5e-7, "floating_gaussian_variational_upper_bound", false,
        NaN, NaN, -0.7608465196784,
        NaN,
        "certificate robustness across N on the three-centre operator (coarsest LG cell; " *
        "no float LG comparison exists at this cell -- dense tier only in the sweep)")
    # ---------------- DECISIVE TEST (correctness, not sharpness) ----------------
    # D1: the most direct test of the open H3^2+ contradiction. At this cell the float
    # LG lower bound is -0.93640124532559 while an independently verified 120-digit
    # variational UPPER bound is -0.936444734 (best estimate -0.936459577 +- 3.8e-07,
    # four basis families agreeing to 7e-07, integrals cross-checked to 5.8e-14).
    # A rigorous lower bound must lie BELOW any variational upper bound. If the
    # certified lower bound also lands above -0.936444734, that is a correctness
    # failure in the chain, not a width problem. Nothing is tuned either way.
    tag == "D1" && return Config("D1","H3plus",2.8,12.0,64,80,32,48,32,
        -0.936459577, 3.8e-7,
        "corrected extended-precision floating-Gaussian best estimate; companion variational upper bound -0.936444734",
        false,
        9.6964832039359e-05, -0.9364012453255908, -0.9363042804935516,
        -0.936444734,
        "DECISIVE TEST: does the CERTIFIED lower bound respect an independently verified variational upper bound?")
    # D2: the smallest-float-margin d=4 cell (+4.62e-06 to the provisional reference).
    # Tests whether the same failure appears at the geometry certified for the paper.
    tag == "D2" && return Config("D2","H3plus",4.0,20.0,64,80,32,48,32,
        -0.76207999734, 3.8e-7,
        "PROVISIONAL floating-Gaussian value; error bar widened from 5e-08 to 3.8e-07 in the Phase-0 correction",
        false,
        0.0003742751312778, -0.7620846156614589, -0.7617103405301809,
        NaN,
        "DECISIVE TEST: smallest-float-margin d=4 cell, at the geometry certified for the paper")
    # --------- sigma_loc rerun: Omega_1, the box the shift parameters are specified for -------
    # Stage A at N_stageA = 80 (the separation certificate L_2 > U_1 fails at 64 with the new
    # shift).  The mu2 separator stage stays at 32/48: raising it to match Stage A would need an
    # oee auxiliary at N' >= 96, dim 115248, ~212 GB -- not available on this host.
    tag == "SL1" && return Config("SL1","H2plus",2.0,10.0,64,80,32,48,80,
        -0.5513170, 0.0, "exact_prolate_spheroidal_half_of_-1.1026342", false,
        NaN, NaN, NaN,
        NaN,
        "sigma_loc rerun on Omega_1: rigorous closed-form coercivity shift, Stage A at N=80")
    tag == "C3b" && error("CONFIG C3b (H3^2+ N=48/N'=64) was CANCELLED by the second scope " *
                          "change; the deliverable is C3a (N=32/N'=48) and the robustness " *
                          "cell is C3c (N=24/N'=40).")
    tag == "C3" && error("CONFIG C3 (H3^2+ N=64/N'=80) was CANCELLED by a scope change; " *
                         "use C3a (N=32/N'=48 de-risk) then C3b (N=48/N'=64 deliverable).")
    tag == "SELFTEST" && return Config("SELFTEST","H2plus",2.0,10.0,32,48,32,48,32,
        -0.5513170, 0.0, "exact_prolate_spheroidal_half_of_-1.1026342", true,
        NaN, NaN, NaN,
        NaN,
        "reproduce the published Omega1 N=32/N'=48 certified numbers")
    error("unknown CONFIG $tag")
end

centres_of(c::Config) = c.system == "H2plus" ? CC.h2plus_centres(2.0) : CC.h3plus_centres(c.d)

iv2(x::IV) = [inf(x), sup(x)]

# --- memory profiling -------------------------------------------------------------------
function meminfo()
    rss = hwm = NaN
    try
        for ln in eachline("/proc/self/status")
            startswith(ln, "VmRSS:") && (rss = parse(Float64, split(ln)[2])/1048576.0)
            startswith(ln, "VmHWM:") && (hwm = parse(Float64, split(ln)[2])/1048576.0)
        end
    catch; end
    return rss, hwm
end
const MEMPROF = Ref(Vector{Any}())
function memlog!(label::String)
    rss, hwm = meminfo()
    prev = isempty(MEMPROF[]) ? 0.0 : MEMPROF[][end]["hwm_gb"]
    push!(MEMPROF[], Dict{String,Any}("step"=>label,"rss_gb"=>rss,"hwm_gb"=>hwm,
                                      "hwm_delta_gb"=>hwm-prev))
    @printf("  [mem] %-34s live=%7.2f GB  peak=%7.2f GB  (peak +%.2f)\n", label, rss, hwm, hwm-prev)
    flush(stdout)
    return nothing
end

# JSON has no NaN/Inf literal: a non-finite Float64 anywhere in the output dict makes
# JSON.print fail and leaves a 0-byte file.  Configs legitimately carry NaN (e.g. C3c
# has no float LG comparison), so map every non-finite number to null recursively.
jsan(x) = x
jsan(x::Float64) = isfinite(x) ? x : nothing
jsan(d::Dict) = Dict(k => jsan(v) for (k,v) in d)
jsan(a::AbstractVector) = [jsan(v) for v in a]
rss() = CC.peak_rss_gb()

# rho (below lambda_k) and sigma (between lambda_k and lambda_{k+1}) for the
# Lehmann-Behnke enclosure of index k, derived from the approximate spectrum so
# that no H2+-specific constant is baked in.
function auto_rs(vals::Vector{Float64}, k::Int)
    rho = vals[k] - 0.05*max(1.0, abs(vals[k])) - 0.02
    sig = 0.5*(vals[k] + vals[k+1])
    return rho, sig
end

# Verified Galerkin enclosure of eigenvalue index k of Hi (float mid Hm supplied).
# `Hm` is passed by REFERENCE in a 1-element vector so this function can drop the caller's
# only reference to the float midpoint before calling lehmann_behnke.  That matters: at
# N=80 the midpoint is 35.4 GB and lehmann_behnke internally allocates TWO more full-size
# dense interval matrices (a promote copy, then the shifted `Ai .- lam*Bi` result), so
# holding the midpoint live across the verification adds 35 GB to the peak for nothing.
function galerkin_encl(Hmref::Vector{Matrix{Float64}}, Hi::Matrix{IV}, D::Int, k::Int;
                       nev::Int=4, rs::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                       free_mid::Bool=true)
    vals, vecs = eigs(Symmetric(Hmref[1]); nev=nev, which=:SR, maxiter=3000)
    memlog!("after eigs (float midpoint live)")
    if free_mid
        Hmref[1] = Matrix{Float64}(undef,0,0); GC.gc()
        memlog!("float midpoint freed")
    end
    p = sortperm(real(vals)); vals = real(vals[p]); vecs = real(vecs[:,p])
    rho, sig = rs === nothing ? auto_rs(vals, k) : rs
    V = reshape(vecs[:,k], D, 1)
    Bi = interval.(sparse(1.0I, D, D))
    t = @elapsed lg = CC.lehmann_behnke(Hi, Bi, vals, V, rho, sig, 1.0, k, k; do_shift=true)
    memlog!("after lehmann_behnke verification")
    return lg[1], vals, vecs, rho, sig, t
end

# =============================================================== STAGE A ======
function stage_A(c::Config, geo, nodes, wts, prefix)
    T0 = time()
    cen, LX, LY, LZ, npan, ts = geo
    # CERT_NSTAGEA overrides the configured Stage-A resolution.  Needed because N_stageA=80
    # OOM-kills on this host: the eee block is D=68921 (76 GB interval) and stage_A still
    # materialises a 38 GB float midpoint for the test-vector eigensolve, so it peaks at
    # 167 GB against ~158 GB available.  The LG stage was restructured in an earlier phase to
    # compute the midpoint during the matvec; stage_A was not.
    NA = haskey(ENV,"CERT_NSTAGEA") ? parse(Int, ENV["CERT_NSTAGEA"]) : c.N_stageA
    rho_box = SigmaLoc.box_rho(cen, LX, LY, LZ)
    Ztot    = SigmaLoc.total_charge(cen)
    @printf("\n===== [%s] STAGE A  Neumann Galerkin enclosures + sharp-L, N=%d =====\n", c.label, NA)
    @printf("  [geometry] box=[-%g,%g]x[-%g,%g]x[-%g,%g]  rho=min dist(nucleus,boundary)=%g  Ztot=%g\n",
            LX,LX,LY,LY,LZ,LZ, rho_box, Ztot)
    flush(stdout)
    t1 = @elapsed (P, K, Deee) = CC.assemble_PK(NA,:even,:even,:even; centres=cen,
            LX=LX, LY=LY, LZ=LZ, npanel=npan, t_star=ts, nodes=nodes, wts=wts)
    @printf("  [eee assemble] D=%d (%.1f s, peak %.1f GB)\n", Deee, t1, rss()); flush(stdout)
    memlog!("eee assembly done (interval matrix live)")

    sv = CC.form_H_inplace!(P, K)
    memlog!("eee form_H_inplace!")
    Hmv = [CC.mid_matrix(P)]
    memlog!("eee mid_matrix formed")
    rs = c.published_rs ? (-0.60,-0.30) : nothing
    mu1, v1, _, r1, s1, tr1 = galerkin_encl(Hmv, P, Deee, 1; rs=rs)
    Hmv = nothing; GC.gc()
    memlog!("eee mu1 enclosure complete")
    CC.restore_diag!(P, sv)
    @printf("  [eee mu1^N] [%.14f,%.14f] w=%.3e  (rho=%.4f sigma=%.4f, %.1f s)\n",
            inf(mu1), sup(mu1), diam(mu1), r1, s1, tr1); flush(stdout)

    # --- the SHIFT ---------------------------------------------------------------------------
    # Two admissible sources, selected by CERT_CEPS_FILE:
    #   (a) unset  -> the analytic closed form sigma_loc(eps)          [v1, pessimistic]
    #   (b) set    -> the certified eta-based Ceps_cert from run_ceps.jl [v2, sharper]
    # Both are UPPER bounds on Ceps^opt = -eta.  The v0 shift -eta_{1,N} is NOT, and is not
    # reachable from here.  When (b) is used the value is read as a POINT interval at the
    # certified upper bound, so no value below Ceps^opt can be consumed.
    Ceps = Dict{Float64,IV}(); slrec = Dict{String,Any}()
    cepsfile = get(ENV, "CERT_CEPS_FILE", "")
    shift_source = cepsfile == "" ? "analytic_sigma_loc" : "certified_eta_based"
    if cepsfile != ""
        cj = JSON.parsefile(cepsfile)
        for e in EPS_GRID
            k = string(e)
            haskey(cj["results"], k) || error("CERT_CEPS_FILE $cepsfile has no eps=$k entry")
            rr = cj["results"][k]
            Ceps[e] = interval(Float64(rr["Ceps_cert"]))
            slrec[k] = Dict{String,Any}("source"=>"certified_eta_based",
                "Ceps_cert"=>rr["Ceps_cert"], "alpha"=>rr["alpha"], "delta"=>rr["delta"],
                "eps_star"=>rr["eps_star"], "C_star"=>rr["C_star"],
                "nu1N_enclosure"=>rr["nu1N_enclosure"], "eta_1N_certified"=>rr["eta_1N_certified"],
                "nu1_lower_certified"=>rr["nu1_lower_certified"],
                "N_aux"=>rr["N_aux"], "file"=>cepsfile,
                "analytic_comparator"=>rr["sigma_loc_analytic_comparator"],
                "improvement_over_analytic"=>rr["improvement_over_analytic"],
                "admissibility"=>rr["admissibility"])
            @printf("  [Ceps certified eps=%.2f] %.14f  (N_aux=%d, alpha=%.6f) vs analytic %.6f -> %.4gx smaller\n",
                    e, sup(Ceps[e]), Int(rr["N_aux"]), Float64(rr["alpha"]),
                    Float64(rr["sigma_loc_analytic_comparator"]),
                    Float64(rr["improvement_over_analytic"])); flush(stdout)
        end
    else
        for e in EPS_GRID
            dl = SIGMA_LOC_DELTA[e]
            r1 = rho_box/2; r2 = rho_box
            Ceps[e] = SigmaLoc.sigma_loc_iv(e, dl, r1, r2, rho_box, Ztot)
            slrec[string(e)] = SigmaLoc.sigma_loc_terms(e, dl, r1, r2, rho_box, Ztot)
            @printf("  [sigma_loc eps=%.2f] delta=%s r1=%g r2=%g rho=%g Ztot=%g -> [%.16f,%.16f] w=%.3e\n",
                    e, string(dl), r1, r2, rho_box, Ztot, inf(Ceps[e]), sup(Ceps[e]), diam(Ceps[e]))
            flush(stdout)
        end
    end

    # --- the auxiliary Ritz value, kept as a DIAGNOSTIC ONLY (never a shift) ---------------
    # Computed as a float Ritz value, not a verified enclosure: it is explicitly not a bound,
    # so paying for lehmann_behnke on a D=68921 matrix twice would buy nothing.
    etarec = Dict{String,Float64}(); t_eta = 0.0
    for e in EPS_GRID
        sv2 = CC.form_H_inplace!(P, K; kin_scale=interval(e))
        Hm2 = CC.mid_matrix(P)
        memlog!(@sprintf("eta diag eps=%.2f mid formed", e))
        te = @elapsed begin
            vals, _ = eigs(Symmetric(Hm2); nev=1, which=:SR, maxiter=3000)
            etarec[string(e)] = minimum(real(vals))
        end
        Hm2 = nothing; GC.gc()
        memlog!(@sprintf("eta diag eps=%.2f done", e))
        CC.restore_diag!(P, sv2)
        t_eta += te
        @printf("  [eta DIAGNOSTIC eps=%.2f] eta_1N_ritz=%.14f  =>  -eta=%.14f\n",
                e, etarec[string(e)], -etarec[string(e)])
        @printf("      NOT a valid shift (Rayleigh-Ritz bounds eta from ABOVE, so -eta_1N <= C_eps^opt); sigma_loc=%.14f used instead\n",
                inf(Ceps[e])); flush(stdout)
    end
    P = nothing; K = nothing; GC.gc()
    memlog!("eee block released")

    t2 = @elapsed (Po, Ko, Doee) = CC.assemble_PK(NA,:odd,:even,:even; centres=cen,
            LX=LX, LY=LY, LZ=LZ, npanel=npan, t_star=ts, nodes=nodes, wts=wts)
    @printf("  [oee assemble] D=%d (%.1f s, peak %.1f GB)\n", Doee, t2, rss()); flush(stdout)
    memlog!("oee assembly done")
    CC.form_H_inplace!(Po, Ko)
    Hmov = [CC.mid_matrix(Po)]
    memlog!("oee mid_matrix formed")
    rso = c.published_rs ? (-0.40,-0.25) : nothing
    mu2, _, _, r2, s2, tr2 = galerkin_encl(Hmov, Po, Doee, 1; rs=rso)
    Hmov = nothing; Po = nothing; Ko = nothing; GC.gc()
    memlog!("oee mu2 enclosure complete, block released")
    @printf("  [oee mu2^N] [%.14f,%.14f] w=%.3e  (rho=%.4f sigma=%.4f)\n",
            inf(mu2), sup(mu2), diam(mu2), r2, s2); flush(stdout)

    bL1 = nothing; be1 = NaN; bL2 = nothing; be2 = NaN
    for e in EPS_GRID
        L = CC.sharp_L_iv(mu1, mu1, NA, e, Ceps[e], LX); L === nothing && continue
        (bL1 === nothing || inf(L) > inf(bL1)) && (bL1 = L; be1 = e)
    end
    for e in EPS_GRID
        L = CC.sharp_L_iv(mu1, mu2, NA, e, Ceps[e], LX); L === nothing && continue
        (bL2 === nothing || inf(L) > inf(bL2)) && (bL2 = L; be2 = e)
    end
    @printf("  [sharp-L] L1=[%.14f,%.14f] w=%.3e (eps=%.2f)\n", inf(bL1), sup(bL1), diam(bL1), be1)
    @printf("  [sharp-L] L2=[%.14f,%.14f] w=%.3e (eps=%.2f)\n", inf(bL2), sup(bL2), diam(bL2), be2)
    flush(stdout)

    # --- separation certificate L_2 > U_1, the test that forced N_stageA = 80 -------------
    U1 = sup(mu1)
    sep_margin = inf(bL2) - U1
    sep_ok = sep_margin > 0
    @printf("  [separation] inf(L_2)=%.14f  vs  U_1=sup(mu_1N)=%.14f   margin=%+.6e   L_2 > U_1: %s\n",
            inf(bL2), U1, sep_margin, sep_ok)
    # Goerisch shift admissibility, verified rather than assumed
    c_prov = max(1.0, ceil((-inf(bL1)+0.25)*4)/4)
    @printf("  [shift check] inf(L_1)=%.14f -> c=1 gives lambda_min(H+I) >= %+.14f (%s); adaptive rule would give c=%.2f\n",
            inf(bL1), 1.0+inf(bL1), 1.0+inf(bL1) > 0 ? "ADMISSIBLE" : "INADMISSIBLE", c_prov)
    flush(stdout)

    rec = Dict{String,Any}("stage"=>"A","N"=>NA,"N_stageA"=>NA,"N_stageA_configured"=>c.N_stageA,"N_stageA_overridden"=>(NA != c.N_stageA),"Deee"=>Deee,"Doee"=>Doee,
        "mu1N_galerkin"=>iv2(mu1),"mu1N_width"=>diam(mu1),
        "mu2N_galerkin"=>iv2(mu2),"mu2N_width"=>diam(mu2),
        "shift_description"=>(shift_source == "certified_eta_based" ?
            ("CERTIFIED eta-based shift: Ceps_cert = -eps*inf(nu_1 lower bound), where nu_1 = " *
             "lambda_min(B) with B = A_eps/eps = -Lap + V/eps. ADMISSIBLE because the projection " *
             "bound is applied to an operator whose Laplacian coefficient is normalised to 1. " *
             "Distinct from the INADMISSIBLE raw -eta_1N: Rayleigh-Ritz bounds eta from above on " *
             "the cosine subspace of H^1(Omega), so -eta_1N <= C_eps^opt at every N. Here the " *
             "verified ENCLOSURE of eta_1N feeds the projection bound rather than being used as " *
             "the shift directly.") :
            ("RIGOROUS closed-form sigma_loc (manuscript lemma) -- replaces the former auxiliary " *
             "Rayleigh-Ritz value -eta_1N, which is inadmissible because Rayleigh-Ritz bounds eta " *
             "from above, so -eta_1N <= C_eps^opt at every N.")),
        "sigma_loc"=>slrec,
        "geometry"=>Dict("LX"=>LX,"LY"=>LY,"LZ"=>LZ,"rho"=>rho_box,"Ztot"=>Ztot,
            "rho_rule"=>"rho = min over nuclei of dist(nucleus, boundary), computed from the geometry record",
            "r1_r2_rule"=>"r1 = rho/2, r2 = rho"),
        "eta_1N_ritz_diagnostic_not_a_bound"=>Dict{String,Any}(
            "values"=>etarec,
            "COMMENT"=>"DIAGNOSTIC ONLY -- NOT A BOUND AND NOT THE SHIFT. eta_1N is the auxiliary " *
                "Rayleigh-Ritz value lambda_1(eps*K + P) on the cosine space, which is a subspace of " *
                "H^1(Omega); Rayleigh-Ritz therefore bounds it from ABOVE, eta_1N >= eta, so " *
                "-eta_1N <= C_eps^opt and using it DIRECTLY as a coercivity shift is INADMISSIBLE. " *
                "Reported in the paper only as an estimate that locates the optimum. The shift " *
                "actually used is under the 'sigma_loc' key with its source named by " *
                "'shift_source'; when that is 'certified_eta_based' the shift derives from a " *
                "VERIFIED ENCLOSURE of eta_1N pushed through the projection bound on the " *
                "normalised operator B = A_eps/eps, which is admissible -- not from this float " *
                "Ritz value.",
            "computed_as"=>"float Ritz value from eigs(Symmetric(mid(H_aux)); nev=1, which=:SR); NOT a " *
                "verified enclosure, since it is explicitly not a bound",
            "minus_eta_would_have_been"=>Dict(k=>-v for (k,v) in etarec),
            "shift_actually_used"=>Dict(string(e)=>sup(Ceps[e]) for e in EPS_GRID),
            "shift_source"=>shift_source),
        "separation_certificate"=>Dict("inf_L2"=>inf(bL2),"U1_sup_mu1N"=>U1,
            "margin"=>sep_margin,"L2_gt_U1"=>sep_ok,
            "note"=>("L_2 > U_1 is required for the separator to separate. Under the analytic " *
                     "sigma_loc shift this failed at N_stageA=64 and forced 80 (72 used in " *
                     "practice); under the certified eta-based shift the margin is far larger. " *
                     "Verified here, not assumed.")),
        "goerisch_shift_admissibility"=>Dict("inf_L1"=>inf(bL1),
            "lambda_min_at_c_1"=>1.0+inf(bL1),"admissible_at_c_1"=>1.0+inf(bL1) > 0,
            "adaptive_rule_would_give"=>c_prov),
        "shift_source"=>shift_source,
        "Ceps"=>Dict(string(e)=>iv2(Ceps[e]) for e in EPS_GRID),
        "Ceps_width"=>Dict(string(e)=>diam(Ceps[e]) for e in EPS_GRID),
        "L1"=>iv2(bL1),"L1_width"=>diam(bL1),"eps1"=>be1,
        "L2"=>iv2(bL2),"L2_width"=>diam(bL2),"eps2"=>be2,
        "lehmann_behnke_rho_sigma"=>Dict("eee"=>[r1,s1],"oee"=>[r2,s2],
            "source"=>c.published_rs ? "published constants" : "auto from approximate spectrum"),
        "memory_profile"=>copy(MEMPROF[]),
        "t_assemble"=>t1+t2,"t_ritz"=>tr1+tr2,"t_eta"=>t_eta,
        "wall_seconds"=>time()-T0,"peak_rss_gb"=>rss())
    open("$(prefix)_stageA.json","w") do f; JSON.print(f, jsan(rec), 2); end
    @printf("  [checkpoint] %s_stageA.json  (%.0f s)\n", prefix, rec["wall_seconds"]); flush(stdout)
    return rec, mu1, mu2, Ceps, bL1, bL2
end

# ======================================================== CERTIFIED mu2 =======
# CERT_NPRIME_SEP overrides the configured oee separator resolution.  Needed because with the
# sigma_loc shift the N'_sep = 48 separator lands BELOW U_1 and does not separate; the measured
# smallest even value that clears U_1 is 52.  Raised explicitly and recorded, never silently.
nprime_sep(c::Config) = haskey(ENV,"CERT_NPRIME_SEP") ? parse(Int, ENV["CERT_NPRIME_SEP"]) : c.Nprime_sep

function stage_mu2(c::Config, geo, nodes, wts, prefix, Ceps, L1, L2)
    NPS = nprime_sep(c)
    T0 = time()
    cen, LX, LY, LZ, npan, ts = geo
    @printf("\n===== [%s] STAGE mu2  certified oee separator, N=%d N'=%d (separator resolution) =====\n",
            c.label, c.N_sep, NPS); flush(stdout)
    ta = @elapsed (P, K, Db) = CC.assemble_PK(NPS,:odd,:even,:even; centres=cen,
            LX=LX, LY=LY, LZ=LZ, npanel=npan, t_star=ts, nodes=nodes, wts=wts)
    @printf("  [oee aux assemble] D=%d (%.1f s, peak %.1f GB)\n", Db, ta, rss()); flush(stdout)

    sv = CC.form_H_inplace!(P, K)
    Hm = CC.mid_matrix(P)
    Hmr = [Hm]; Hm = nothing
    mu2N, vals, vecs, r2, s2, _ = galerkin_encl(Hmr, P, Db, 1; free_mid=false)
    mu3N, _, _, r3, s3, _       = galerkin_encl(Hmr, P, Db, 2; free_mid=false)
    @printf("  [oee Galerkin] mu2^N'=[%.14f,%.14f] mu3^N'=[%.14f,%.14f]\n",
            inf(mu2N), sup(mu2N), inf(mu3N), sup(mu3N)); flush(stdout)

    bsep = nothing; be = NaN
    for e in EPS_GRID
        L = CC.sharp_L_iv(mu2N, mu3N, NPS, e, Ceps[e], LX); L === nothing && continue
        (bsep === nothing || inf(L) > inf(bsep)) && (bsep = L; be = e)
    end
    sep = inf(bsep); sep_valid = sep > sup(mu2N)
    @printf("  [oee-2nd sharp-L] separator=%.14f (eps=%.2f) valid=%s\n", sep, be, sep_valid)
    flush(stdout)

    pidx = Int[]
    nxs = length(CC.sector_idx(c.N_sep,:odd)); nys = length(CC.sector_idx(c.N_sep,:even))
    nxb = length(CC.sector_idx(NPS,:odd)); nyb = length(CC.sector_idx(NPS,:even))
    for cz in 1:nys, cy in 1:nys, cx in 1:nxs
        push!(pidx, (cz-1)*nyb*nxb + (cy-1)*nxb + cx)
    end
    Hsub = Hmr[1][pidx,pidx]; es = eigen(Symmetric(Hsub), 1:1); vs = vec(es.vectors[:,1])
    Hsub = nothing; Hmr[1] = Matrix{Float64}(undef,0,0); GC.gc()
    memlog!("mu2 stage: primary sub-block taken, float midpoint freed")
    v = zeros(Float64, Db); v[pidx] .= vs; v ./= norm(v); vI = interval.(v)
    CC.restore_diag!(P, sv)

    CC.form_H_inplace!(P, K; sigma=interval(CSHIFT()))
    Kdiag = copy(K); K = nothing; GC.gc()
    Hms = CC.mid_matrix(P)
    Hv = CC.imatvec(P, vI); A0 = CC.idot(vI, Hv)
    w = zeros(Float64, Db)
    tcg = @elapsed cg!(w, Symmetric(Hms), v; reltol=1e-12, maxiter=8000)
    cg_res = norm(Hms*w .- v)/norm(v); Hms = nothing; GC.gc()
    wI = interval.(w); Hw = CC.imatvec(P, wI); r = vI .- Hw
    wv = CC.idot(wI, vI); wHw = CC.idot(wI, Hw); r2sq = CC.idot(r, r)
    lam_min_lb = CSHIFT() + inf(L1)
    @assert lam_min_lb > 0 "shifted operator not certified SPD: c + inf(L1) = $lam_min_lb"
    corr = interval(0.0, sup(r2sq)/lam_min_lb)
    A2_gal = interval(2.0)*wv - wHw + corr
    P = nothing; GC.gc()
    memlog!("mu2 stage: interval operator released before the V^2 form")

    # ---- corrected realization, oee separator stage (same defect, same repair) ----------
    Ixo = CC.sector_idx(NPS, :odd); Iyo = CC.sector_idx(NPS, :even)
    nxo = length(Ixo); nyo = length(Iyo)
    # the closures record what they were actually built with, for the G4 comparison below
    prov_o = Dict{String,Any}("nx_calls"=>0, "nyz_calls"=>0)
    momx_o  = (s,t) -> begin
        prov_o["npanel"] = npan; prov_o["t_star"] = ts; prov_o["LX"] = LX
        prov_o["Nprime_x"] = NPS; prov_o["nx_calls"] = prov_o["nx_calls"] + 1
        CC.cos_moment_matrix_1d(NPS, LX, s, t, nodes, wts; npanel=npan, t_star=ts)[Ixo,Ixo]
    end
    momyz_o = (s,t) -> begin
        prov_o["LYZ"] = LY; prov_o["Nprime_yz"] = NPS
        prov_o["nyz_calls"] = prov_o["nyz_calls"] + 1
        CC.cos_moment_matrix_1d(NPS, LY, s, t, nodes, wts; npanel=npan, t_star=ts)[Iyo,Iyo]
    end
    tg, tw = CC.build_t_grid(48)
    grec, gcands = goerisch_block("oee", wI, Hw, r, Kdiag, wv, wHw, r2sq, CSHIFT(), NPS, LX,
                                  nxo, nyo, nyo, momx_o, momyz_o, cen, tg, tw, false, prov_o)
    gate_G4!("mu2_goerisch",
        ["nx*ny*nz == D_aux"                => nxo*nyo*nyo == Db,
         "V2 sector == oee separator sector"=> (nxo == length(CC.sector_idx(NPS,:odd))),
         "V2 Nprime == separator Nprime"    => (grec["Nprime_used"] == NPS &&
                                                get(prov_o,"Nprime_x",-1) == NPS &&
                                                get(prov_o,"Nprime_yz",-1) == NPS),
         "ceps geometry: LX"                => (isempty(GCEPSREC) || GCEPSREC["_box"]["LX"] == LX),
         "ceps geometry: Ztot"              => (isempty(GCEPSREC) || GCEPSREC["_box"]["Ztot"] == sum(x[2] for x in cen)),
         "ceps geometry: system"            => (isempty(GCEPSREC) || GCEPSREC["_system"] == c.system),
         "ceps geometry: d_paper"           => (isempty(GCEPSREC) || GCEPSREC["_d_paper"] == c.d),
         "npanel/t_star from this run"      => (grec["npanel_used"] == npan &&
                                                grec["t_star_used"] == ts &&
                                                grec["LX_used"] == LX && grec["LYZ_used"] == LY &&
                                                grec["momx_calls"] > 0 && grec["momyz_calls"] > 0)],
        Dict{String,Any}("D_aux"=>Db,"nx"=>nxo,"ny"=>nyo,"Nprime_sep"=>NPS,"LX"=>LX,"LY"=>LY,
            "npanel"=>npan,"t_star"=>ts,"nt"=>48,"centres"=>string(cen),
            "ceps_file"=>get(GCEPSREC,"_file",""),"ceps_N_aux"=>get(GCEPSREC,"_N_aux",nothing)))
    Hw = nothing; Kdiag = nothing; GC.gc()

    rho_p = interval(sep) + interval(CSHIFT())
    function mobius_mu2(A2x::IV)
        Ax = A0 - rho_p; Bx = A0 - interval(2.0)*rho_p + rho_p*rho_p*A2x; nux = Ax/Bx
        Lx = (rho_p - rho_p/(interval(1.0)-nux)) - interval(CSHIFT())
        return Lx, Bx, nux
    end
    gvar = Dict{String,Any}(); bestg = nothing
    for (e, ce, A2G, cr, df, s2) in gcands
        Lx, Bx, nux = mobius_mu2(A2G)
        okx = inf(Bx) > 0 && sup(nux) < 1
        gvar[string(e)] = Dict{String,Any}("eps_prime"=>e, "Ceps_cert"=>sup(ce),
            "A2_goerisch"=>[inf(A2G),sup(A2G)], "correction_sup"=>sup(cr),
            "B"=>[inf(Bx),sup(Bx)], "nu"=>[inf(nux),sup(nux)],
            "B_positive"=>inf(Bx)>0, "nu_lt_1"=>sup(nux)<1, "certificate_valid"=>okx,
            "L2_LG"=>[isfinite(inf(Lx)) ? inf(Lx) : nothing,
                      isfinite(sup(Lx)) ? sup(Lx) : nothing])
        @printf("  [oee LG corrected] eps'=%.4g -> mu2(Omega) >= %.14f  B>0=%s nu<1=%s\n",
                e, inf(Lx), inf(Bx)>0, sup(nux)<1); flush(stdout)
        okx && (bestg === nothing || inf(Lx) > bestg[2]) && (bestg = (e, inf(Lx), A2G))
    end
    if GOFF || bestg === nothing
        A2 = A2_gal; eps_used = nothing
        GOFF || println("  [oee LG corrected] NO ADMISSIBLE eps' -- the uncorrected A_2 is used " *
                        "below FOR DIAGNOSIS ONLY; the separator it yields is NOT certified")
    else
        eps_used = bestg[1]; A2 = bestg[3]
    end
    L2_LG, B, nu = mobius_mu2(A2)
    A = A0 - rho_p
    Bpos = inf(B) > 0; nult1 = sup(nu) < 1
    @printf("  [oee LG] A0=[%.14f,%.14f] A2=[%.14f,%.14f] cg_res=%.2e\n",
            inf(A0), sup(A0), inf(A2), sup(A2), cg_res)
    @printf("  [oee LG] B=[%.4e,%.4e] pos=%s  nu=[%.5f,%.5f] <1=%s\n",
            inf(B), sup(B), Bpos, inf(nu), sup(nu), nult1)
    @printf("  [oee LG] mu2(Omega) >= %.14f   (width %.3e)\n", inf(L2_LG), diam(L2_LG))
    flush(stdout)

    fin(x) = isfinite(x) ? x : nothing
    rec = Dict{String,Any}("stage"=>"mu2","N"=>c.N_sep,"Nprime"=>NPS,"Nprime_sep_configured"=>c.Nprime_sep,"Nprime_sep_overridden"=>(NPS != c.Nprime_sep),"D_aux"=>Db,"c"=>CSHIFT(),
        "mu2N_galerkin"=>iv2(mu2N),"mu2N_width"=>diam(mu2N),
        "mu3N_galerkin"=>iv2(mu3N),"mu3N_width"=>diam(mu3N),
        "oee2nd_sharpL_separator"=>sep,"sep_eps"=>be,"sep_valid"=>sep_valid,
        "A0"=>iv2(A0),"A0_width"=>diam(A0),"A2"=>iv2(A2),"A2_width"=>diam(A2),
        "A2_residual_corr_hi"=>sup(corr),"r2_hi"=>sup(r2sq),"lam_min_lb"=>lam_min_lb,
        "A2_galerkin"=>iv2(A2_gal),"A2_galerkin_width"=>diam(A2_gal),
        "goerisch"=>grec,"goerisch_variants"=>gvar,"eps_prime_used"=>eps_used,
        "corrected_realization_used"=>(eps_used !== nothing && !GOFF),
        "B"=>[fin(inf(B)),fin(sup(B))],"nu"=>[fin(inf(nu)),fin(sup(nu))],
        "B_positive"=>Bpos,"nu_lt_1"=>nult1,"certificate_valid"=>(Bpos && nult1),
        "cg_residual"=>cg_res,"t_cg"=>tcg,
        "mu2_Omega_lower_certified"=>fin(inf(L2_LG)),
        "L2_LG"=>[fin(inf(L2_LG)),fin(sup(L2_LG))],"L2_LG_width"=>fin(diam(L2_LG)),
        "lehmann_behnke_rho_sigma"=>Dict("mu2"=>[r2,s2],"mu3"=>[r3,s3]),
        "t_assemble"=>ta,"wall_seconds"=>time()-T0,"peak_rss_gb"=>rss())
    open("$(prefix)_stagemu2.json","w") do f; JSON.print(f, jsan(rec), 2); end
    @printf("  [checkpoint] %s_stagemu2.json  (%.0f s)\n", prefix, rec["wall_seconds"]); flush(stdout)
    return rec, (Bpos && nult1) ? inf(L2_LG) : nothing
end

# ============================================== eee LEHMANN-GOERISCH LOWER ====
function stage_lg(c::Config, geo, nodes, wts, prefix, L1, L2, mu2_cert)
    T0 = time()
    cen, LX, LY, LZ, npan, ts = geo
    # best available certified separator: the certified mu2(Omega) lower bound if
    # its certificate held, else the Stage-A L2.  Both are rigorous lower bounds
    # on lambda_2; the larger one is the sharper separator.
    cands = Tuple{String,Float64}[("stageA_L2", inf(L2))]
    mu2_cert !== nothing && push!(cands, ("certified_oee_mu2", mu2_cert))
    src, rho = cands[argmax([x[2] for x in cands])]
    @printf("\n===== [%s] STAGE LG  eee lower bound, N=%d N'=%d, rho=%.14f (%s) =====\n",
            c.label, c.N, c.Nprime, rho, src); flush(stdout)

    ta = @elapsed (P, K, Db) = CC.assemble_PK(c.Nprime,:even,:even,:even; centres=cen,
            LX=LX, LY=LY, LZ=LZ, npanel=npan, t_star=ts, nodes=nodes, wts=wts)
    @printf("  [eee aux assemble] D=%d (%.1f s, peak %.1f GB)\n", Db, ta, rss()); flush(stdout)

    CC.form_H_inplace!(P, K; sigma=interval(CSHIFT()))
    Kdiag = copy(K); K = nothing; GC.gc(true)

    nxs = length(CC.sector_idx(c.N,:even)); nxb = length(CC.sector_idx(c.Nprime,:even))
    pidx = Int[]
    for cz in 1:nxs, cy in 1:nxs, cx in 1:nxs
        push!(pidx, (cz-1)*nxb*nxb + (cy-1)*nxb + cx)
    end
    # (D7): slice the primary block straight out of the interval operator and run
    # CG against a midpoint-on-the-fly operator, so the 38 GB float midpoint is
    # never allocated.
    Hsub = CC.mid_subblock(P, pidx)
    es = eigen(Symmetric(Hsub), 1:1)
    mu1_primary_shifted = es.values[1]; vs = vec(es.vectors[:,1])
    Hsub = nothing; es = nothing; GC.gc(true)
    @printf("  [primary block] %d^3 = %d inside aux %d^3 = %d, peak %.1f GB\n",
            nxs, length(pidx), nxb, Db, rss()); flush(stdout)
    v = zeros(Float64, Db); v[pidx] .= vs; v ./= norm(v); vI = interval.(v)

    Hv = CC.imatvec(P, vI); A0 = CC.idot(vI, Hv)
    Aop = CC.MidOp(P)
    tcg = @elapsed ((w, nit, cg_res) = CC.cg_mid(Aop, v; reltol=1e-12, maxiter=8000))
    @printf("  [cg] %d iters, rel residual %.2e (%.1f s)\n", nit, cg_res, tcg); flush(stdout)
    wI = interval.(w); Hw = CC.imatvec(P, wI); r = vI .- Hw
    wv = CC.idot(wI, vI); wHw = CC.idot(wI, Hw); r2sq = CC.idot(r, r)
    P = nothing; GC.gc()
    memlog!("lg stage: interval operator released before the V^2 form")

    lam_min_lb = CSHIFT() + inf(L1)
    @assert lam_min_lb > 0 "shifted operator not certified SPD: c + inf(L1) = $lam_min_lb"
    corr = interval(0.0, sup(r2sq)/lam_min_lb)
    A2_gal = interval(2.0)*wv - wHw + corr

    # ---- corrected realization, eee primary stage ---------------------------------------
    Ie = CC.sector_idx(c.Nprime, :even); nie = length(Ie)
    prov_e = Dict{String,Any}("nx_calls"=>0, "nyz_calls"=>0)
    momx_e  = (s,t) -> begin
        prov_e["npanel"] = npan; prov_e["t_star"] = ts; prov_e["LX"] = LX
        prov_e["Nprime_x"] = c.Nprime; prov_e["nx_calls"] = prov_e["nx_calls"] + 1
        CC.cos_moment_matrix_1d(c.Nprime, LX, s, t, nodes, wts; npanel=npan, t_star=ts)[Ie,Ie]
    end
    momyz_e = (s,t) -> begin
        prov_e["LYZ"] = LY; prov_e["Nprime_yz"] = c.Nprime
        prov_e["nyz_calls"] = prov_e["nyz_calls"] + 1
        CC.cos_moment_matrix_1d(c.Nprime, LY, s, t, nodes, wts; npanel=npan, t_star=ts)[Ie,Ie]
    end
    tg, tw = CC.build_t_grid(48)
    grec, gcands = goerisch_block("eee", wI, Hw, r, Kdiag, wv, wHw, r2sq, CSHIFT(), c.Nprime,
                                  LX, nie, nie, nie, momx_e, momyz_e, cen, tg, tw, true, prov_e)
    gate_G4!("lg_goerisch",
        ["nx*ny*nz == D_aux"            => nie*nie*nie == Db,
         "V2 Nprime == LG Nprime"       => (grec["Nprime_used"] == c.Nprime &&
                                            get(prov_e,"Nprime_x",-1) == c.Nprime &&
                                            get(prov_e,"Nprime_yz",-1) == c.Nprime),
         "ceps geometry: LX"            => (isempty(GCEPSREC) || GCEPSREC["_box"]["LX"] == LX),
         "ceps geometry: Ztot"          => (isempty(GCEPSREC) || GCEPSREC["_box"]["Ztot"] == sum(x[2] for x in cen)),
         "ceps geometry: system"        => (isempty(GCEPSREC) || GCEPSREC["_system"] == c.system),
         "ceps geometry: d_paper"       => (isempty(GCEPSREC) || GCEPSREC["_d_paper"] == c.d),
         "npanel/t_star from this run"  => (grec["npanel_used"] == npan &&
                                            grec["t_star_used"] == ts &&
                                            grec["LX_used"] == LX && grec["LYZ_used"] == LY &&
                                            grec["momx_calls"] > 0 && grec["momyz_calls"] > 0)],
        Dict{String,Any}("D_aux"=>Db,"n_per_axis"=>nie,"Nprime"=>c.Nprime,"LX"=>LX,"LY"=>LY,
            "npanel"=>npan,"t_star"=>ts,"nt"=>48,"centres"=>string(cen),
            "ceps_file"=>get(GCEPSREC,"_file",""),"ceps_N_aux"=>get(GCEPSREC,"_N_aux",nothing)))
    Hw = nothing; Kdiag = nothing; GC.gc()

    # Hypothesis of Theorem lg-sharp: rho must exceed Lambda_n, the Rayleigh quotient of the
    # test vector.  A1 = <v,v> = 1 by construction, so Lambda_n = A0 - c.  Checked BEFORE the
    # pencil is formed, and recorded either way.
    Lambda_n = sup(A0) - CSHIFT()
    rho_gt_Lambda_n = rho > Lambda_n
    @printf("  [hypothesis] Lambda_n = sup(A0) - c = %.14f ; rho = %.14f ; rho > Lambda_n: %s (margin %+.6e)\n",
            Lambda_n, rho, rho_gt_Lambda_n, rho - Lambda_n); flush(stdout)
    @assert rho_gt_Lambda_n "Theorem lg-sharp hypothesis violated: rho=$rho <= Lambda_n=$Lambda_n"

    # A0 and A2 do not depend on rho, so the Moebius step is evaluated at EVERY
    # candidate separator from the one assembly.  This is what lets a single file
    # carry both the sharpest bound and the published-separator variant.
    function mobius(rr::Float64, A2x::IV)
        rp = interval(rr) + interval(CSHIFT())
        Ax = A0 - rp
        Bx = A0 - interval(2.0)*rp + rp*rp*A2x
        nux = Ax/Bx
        Lx = (rp - rp/(interval(1.0)-nux)) - interval(CSHIFT())
        return Lx, Bx, nux
    end
    # Every admissible eps' yields a VALID certificate; take the sharpest, exactly as this
    # driver already takes the sharpest separator and the sharpest sharp-L eps.  No
    # parameter is being tuned to make an inequality hold: inadmissible eps' are discarded
    # by the c > Ceps' test, not by their effect on the bound.
    gvar = Dict{String,Any}(); bestg = nothing
    for (e, ce, A2G, cr, df, s2) in gcands
        Lx, Bx, nux = mobius(rho, A2G)
        okx = inf(Bx) > 0 && sup(nux) < 1
        gvar[string(e)] = Dict{String,Any}("eps_prime"=>e,"Ceps_cert"=>sup(ce),
            "A2_goerisch"=>iv2(A2G),"correction_sup"=>sup(cr),
            "defect_sup"=>sup(df),"sigma2_const_mode_sup"=>sup(s2),
            "B"=>[inf(Bx),sup(Bx)],"nu"=>[inf(nux),sup(nux)],
            "B_positive"=>inf(Bx)>0,"nu_lt_1"=>sup(nux)<1,"certificate_valid"=>okx,
            "L1_LG"=>[isfinite(inf(Lx)) ? inf(Lx) : nothing,
                      isfinite(sup(Lx)) ? sup(Lx) : nothing],
            "L1_LG_width"=>isfinite(diam(Lx)) ? diam(Lx) : nothing)
        @printf("  [eee LG corrected] eps'=%.4g Ceps'=%.8f -> lambda_1 >= %.14f  B>0=%s nu<1=%s\n",
                e, sup(ce), inf(Lx), inf(Bx)>0, sup(nux)<1); flush(stdout)
        okx && (bestg === nothing || inf(Lx) > bestg[2]) && (bestg = (e, inf(Lx), A2G))
    end
    if GOFF || bestg === nothing
        A2 = A2_gal; eps_used = nothing
        GOFF || println("  [eee LG corrected] NO ADMISSIBLE eps' -- the uncorrected A_2 is used " *
                        "below FOR DIAGNOSIS ONLY; that bound is NOT certified for the operator")
    else
        eps_used = bestg[1]; A2 = bestg[3]
    end
    variants = Dict{String,Any}()
    for (nm, rr) in cands
        Lx, Bx, nux = mobius(rr, A2)
        variants[nm] = Dict{String,Any}("rho"=>rr,
            "L1_LG"=>[isfinite(inf(Lx)) ? inf(Lx) : nothing, isfinite(sup(Lx)) ? sup(Lx) : nothing],
            "L1_LG_width"=>isfinite(diam(Lx)) ? diam(Lx) : nothing,
            "B"=>[inf(Bx),sup(Bx)],"nu"=>[inf(nux),sup(nux)],
            "B_positive"=>inf(Bx)>0,"nu_lt_1"=>sup(nux)<1,
            "certificate_valid"=>(inf(Bx)>0 && sup(nux)<1))
        @printf("  [eee LG variant %-20s] rho=%.14f -> lower=%.14f  B>0=%s nu<1=%s\n",
                nm, rr, inf(Lx), inf(Bx)>0, sup(nux)<1); flush(stdout)
    end
    L1_LG, B, nu = mobius(rho, A2)
    Bpos = inf(B) > 0; nult1 = sup(nu) < 1
    @printf("  [eee LG] A0=[%.14f,%.14f] A2=[%.14f,%.14f] cg_res=%.2e\n",
            inf(A0), sup(A0), inf(A2), sup(A2), cg_res)
    @printf("  [eee LG] B=[%.4e,%.4e] pos=%s  nu=[%.5f,%.5f] <1=%s\n",
            inf(B), sup(B), Bpos, inf(nu), sup(nu), nult1)
    @printf("  [eee LG] lambda_1 >= %.14f   (L1_LG width %.3e)\n", inf(L1_LG), diam(L1_LG))
    flush(stdout)

    fin(x) = isfinite(x) ? x : nothing
    rec = Dict{String,Any}("stage"=>"lg","N"=>c.N,"Nprime"=>c.Nprime,"D_aux"=>Db,"c"=>CSHIFT(),
        "hypothesis_rho_gt_Lambda_n"=>Dict("Lambda_n"=>Lambda_n,"rho"=>rho,
            "holds"=>rho_gt_Lambda_n,"margin"=>rho-Lambda_n,
            "note"=>"Lambda_n = sup(A0) - c, using A1 = <v,v> = 1; checked before the pencil was formed"),
        "separator_rho"=>rho,"separator_source"=>src,
        "separator_candidates"=>Dict(k=>v for (k,v) in cands),
        "separator_variants"=>variants,
        "mu1_primary_float_unshifted"=>mu1_primary_shifted - CSHIFT(),
        "A0"=>iv2(A0),"A0_width"=>diam(A0),"A2"=>iv2(A2),"A2_width"=>diam(A2),
        "A2_residual_corr_hi"=>sup(corr),"r2_hi"=>sup(r2sq),"lam_min_lb"=>lam_min_lb,
        "A2_galerkin"=>iv2(A2_gal),"A2_galerkin_width"=>diam(A2_gal),
        "goerisch"=>grec,"goerisch_variants"=>gvar,"eps_prime_used"=>eps_used,
        "corrected_realization_used"=>(eps_used !== nothing && !GOFF),
        "B"=>[fin(inf(B)),fin(sup(B))],"nu"=>[fin(inf(nu)),fin(sup(nu))],
        "B_positive"=>Bpos,"nu_lt_1"=>nult1,"certificate_valid"=>(Bpos && nult1),
        "cg_residual"=>cg_res,"t_cg"=>tcg,
        "lambda1_lower_certified"=>fin(inf(L1_LG)),
        "L1_LG"=>[fin(inf(L1_LG)),fin(sup(L1_LG))],"L1_LG_width"=>fin(diam(L1_LG)),
        "t_assemble"=>ta,"wall_seconds"=>time()-T0,"peak_rss_gb"=>rss())
    open("$(prefix)_stagelg.json","w") do f; JSON.print(f, jsan(rec), 2); end
    @printf("  [checkpoint] %s_stagelg.json  (%.0f s)\n", prefix, rec["wall_seconds"]); flush(stdout)
    return rec, (Bpos && nult1) ? inf(L1_LG) : nothing
end

# ================================================ DIRICHLET UPPER BOUND ======
function stage_dir(c::Config, geo, nodes, wts, prefix)
    T0 = time()
    cen, LX, LY, LZ, npan, ts = geo
    @printf("\n===== [%s] STAGE DIR  Dirichlet Rayleigh-Ritz upper, N=%d =====\n", c.label, c.N)
    flush(stdout)
    ta = @elapsed (H, D, nx) = CC.assemble_dirichlet(c.N; centres=cen, LX=LX, LY=LY, LZ=LZ,
            npanel=npan, t_star=ts, nodes=nodes, wts=wts)
    @printf("  [dirichlet assemble] D=%d modes/axis=%d (%.1f s, peak %.1f GB)\n",
            D, nx, ta, rss()); flush(stdout)
    Hm = CC.mid_matrix(H)
    es = eigen(Symmetric(Hm), 1:1)
    fl = es.values[1]; v = vec(es.vectors[:,1]); v ./= norm(v)
    Hm = nothing; GC.gc()
    vI = interval.(v)
    Hv = CC.imatvec(H, vI); num = CC.idot(vI, Hv); den = CC.idot(vI, vI)
    H = nothing; GC.gc()
    RQ = num/den; ub = sup(RQ)
    @printf("  [dirichlet] float lambda1^D_N = %.14f\n", fl)
    @printf("  [dirichlet] RQ=[%.14f,%.14f] w=%.3e => certified UPPER = %.14f\n",
            inf(RQ), sup(RQ), diam(RQ), ub); flush(stdout)
    rec = Dict{String,Any}("stage"=>"dirichlet","N"=>c.N,"D"=>D,"modes_per_axis"=>nx,
        "float_lambda1D_N"=>fl,"RQ_interval"=>iv2(RQ),"RQ_width"=>diam(RQ),
        "mass_den"=>iv2(den),"lambda1_upper_certified"=>ub,
        "t_assemble"=>ta,"wall_seconds"=>time()-T0,"peak_rss_gb"=>rss())
    open("$(prefix)_stagedir.json","w") do f; JSON.print(f, jsan(rec), 2); end
    @printf("  [checkpoint] %s_stagedir.json  (%.0f s)\n", prefix, rec["wall_seconds"]); flush(stdout)
    return rec, ub
end

# ===================================================================== main ==
# ---------- Stage A and the Dirichlet stage read from THIS configuration's own files ------
# Neither stage uses Lehmann-Goerisch, so neither is touched by this repair.  Reading them
# back from the recorded per-stage files of the SAME configuration costs nothing in rigour
# and saves re-running them; every field that could differ is asserted (gate G4).
function load_stageA_from(pfx::String, c::Config, LX, LY, LZ)
    f = "$(pfx)_stageA.json"
    J = JSON.parsefile(f)
    g = J["geometry"]
    gate_G4!("stageA_from_file",
        ["file box LX"      => (g["LX"] == LX),
         "file box LY"      => (g["LY"] == LY),
         "file box LZ"      => (g["LZ"] == LZ),
         "file Ztot"        => (g["Ztot"] == sum(x[2] for x in centres_of(c))),
         "file rho"         => (g["rho"] == g["rho"]),
         "file N_stageA == this run's N_stageA" => (J["N_stageA"] == NSTAGEA(c)),
         "file has L1,L2,Ceps" => (haskey(J,"L1") && haskey(J,"L2") && haskey(J,"Ceps"))],
        Dict{String,Any}("file"=>f,"geometry"=>g,"N_stageA"=>J["N_stageA"],
                         "L1"=>J["L1"],"L2"=>J["L2"],"Ceps"=>J["Ceps"]))
    iv(x) = interval(Float64(x[1]), Float64(x[2]))
    Ceps = Dict{Float64,IV}()
    for (k,v) in J["Ceps"]; Ceps[parse(Float64,String(k))] = iv(v); end
    @printf("  [stage A] READ FROM FILE %s : N_stageA=%d L1=[%.14f,%.14f] L2=[%.14f,%.14f]\n",
            f, J["N_stageA"], J["L1"][1], J["L1"][2], J["L2"][1], J["L2"][2]); flush(stdout)
    rec = Dict{String,Any}("stage"=>"A","source"=>"READ FROM $f (unchanged by this repair)",
        "N_stageA"=>J["N_stageA"], "L1"=>J["L1"], "L2"=>J["L2"], "Ceps"=>J["Ceps"],
        "mu1N_galerkin"=>J["mu1N_galerkin"], "mu2N_galerkin"=>J["mu2N_galerkin"],
        "geometry"=>g)
    return rec, iv(J["mu1N_galerkin"]), iv(J["mu2N_galerkin"]), Ceps, iv(J["L1"]), iv(J["L2"])
end
function load_stagedir_from(pfx::String, c::Config)
    f = "$(pfx)_stagedir.json"
    J = JSON.parsefile(f)
    gate_G4!("stagedir_from_file",
        ["file N == config N" => (J["N"] == c.N),
         "file has certified upper" => haskey(J,"lambda1_upper_certified")],
        Dict{String,Any}("file"=>f,"N"=>J["N"],"upper"=>J["lambda1_upper_certified"]))
    @printf("  [stage dir] READ FROM FILE %s : lambda_1 <= %.14f\n",
            f, J["lambda1_upper_certified"]); flush(stdout)
    rec = Dict{String,Any}("stage"=>"dirichlet","source"=>"READ FROM $f (unchanged by this repair)",
        "N"=>J["N"], "RQ_interval"=>J["RQ_interval"],
        "lambda1_upper_certified"=>J["lambda1_upper_certified"])
    return rec, Float64(J["lambda1_upper_certified"])
end

function main()
    tag = ARGS[1]
    CC.THREADED[] = get(ENV, "CERT_THREADED", "1") == "1"
    haskey(ENV,"CERT_BLAS") && BLAS.set_num_threads(parse(Int, ENV["CERT_BLAS"]))
    c = get_config(tag)
    prefix = "cert_$(tag)"
    T0 = time()
    cen = centres_of(c)
    LX = c.L; LY = 0.8*c.L; LZ = 0.8*c.L
    npan = CC.npanel_rule(LX)
    ts, tail = CC.choose_t_star(cen, LX, LY, LZ; tol=1e-13)
    geo = (cen, LX, LY, LZ, npan, ts)
    sig = CC.sigma_conf(cen, LX, LY, LZ)
    clearance = LX - maximum(abs(x[1][1]) for x in cen)

    @printf("CONFIG=%s system=%s d=%g box=[-%g,%g]x[-%g,%g]^2 N=%d N'=%d\n",
            tag, c.system, c.d, LX, LX, LY, LY, c.N, c.Nprime)
    @printf("npanel=%d (rule round(4.8*LX))  t_star=%.2f (tail %.2e <= 1e-13)  clearance=%.2f\n",
            npan, ts, tail, clearance)
    @printf("sigma(Omega)=%.9f  lambda_ref=%.9f  confinement margin=%.6f  holds=%s\n",
            sig, c.lambda_ref, sig - c.lambda_ref, sig > c.lambda_ref)
    @printf("threads=%d BLAS=%d threaded_interval=%s julia=%s\n", Threads.nthreads(), BLAS.get_num_threads(), CC.THREADED[], VERSION)
    flush(stdout)
    @assert sig > c.lambda_ref "confinement hypothesis fails: sigma=$sig <= lambda_1=$(c.lambda_ref)"

    nodes, wts = CC.MomentsVerified.gl_reference(24)
    println("certified GL reference rule built (n_p=24)."); flush(stdout)

    load_gceps!()
    if isempty(GCEPS)
        println("NO certified C_eps' loaded (set CERT_G_CEPSFILE) -- the corrected realization " *
                "cannot be certified in this run.")
        @assert GOFF "refusing to run: CERT_G_CEPSFILE=$(GCEPSFILE) yielded no certified " *
            "C_eps'. A run that cannot certify the corrected realization must be asked for " *
            "explicitly with CERT_G_OFF=1."
    else
        for (k,v) in sort(collect(GCEPS); by=first)
            @printf("certified C_eps' : eps'=%-6.4g Ceps'=%.10f   (from %s, N_aux=%s)\n",
                    k, sup(v), GCEPSFILE, string(get(GCEPSREC,"_N_aux",nothing)))
        end
    end
    flush(stdout)

    recA, mu1, mu2, Ceps, L1, L2 = GFROM == "" ?
        stage_A(c, geo, nodes, wts, prefix) : load_stageA_from(GFROM, c, LX, LY, LZ)
    cs = set_shift_from_L1!(L1)
    @printf("  [Goerisch shift] inf(L1)=%.14f -> c=%.2f, certified lambda_min(H+cI) >= %.14f %s\n",
            inf(L1), cs, cs + inf(L1), cs > 1.0 ? "(RAISED above the H2+ value 1.0)" : "")
    flush(stdout)
    GC.gc(true); recM, mu2_cert = stage_mu2(c, geo, nodes, wts, prefix, Ceps, L1, L2)
    GC.gc(true); recL, lam_lo   = stage_lg(c, geo, nodes, wts, prefix, L1, L2, mu2_cert)
    GC.gc(true); recD, lam_hi   = GFROM == "" ?
        stage_dir(c, geo, nodes, wts, prefix) : load_stagedir_from(GFROM, c)

    ok = lam_lo !== nothing
    W = ok ? lam_hi - lam_lo : nothing
    contains_ref = ok ? (lam_lo <= c.lambda_ref <= lam_hi) : nothing
    @printf("\n===== [%s] CERTIFIED ENCLOSURE =====\n", tag)
    if ok
        @printf("  lambda_1 in [%.14f, %.14f]   WIDTH = %.6e\n", lam_lo, lam_hi, W)
        @printf("  contains lambda_ref %.9f : %s\n", c.lambda_ref, contains_ref)
    else
        println("  NO CERTIFIED LOWER BOUND -- the Lehmann-Goerisch certificate failed.")
    end
    flush(stdout)

    out = Dict{String,Any}(
        "config"=>tag, "system"=>c.system, "d_paper"=>c.d,
        "box"=>Dict("LX"=>LX,"LY"=>LY,"LZ"=>LZ,
                    "label"=>"[-$(LX),$(LX)] x [-$(LY),$(LY)]^2"),
        "N"=>c.N, "Nprime"=>c.Nprime,
        "separator_resolution"=>Dict("N_sep"=>c.N_sep,"Nprime_sep"=>c.Nprime_sep,
            "note"=>"Stage-A (C_eps, L1, L2) and the certified oee mu2 separator are computed at " *
                    "N_sep/Nprime_sep, not at N/Nprime. Both are rigorous lower bounds on their " *
                    "targets at ANY resolution, so this is sound; it is also the published pattern " *
                    "(the Omega2 N=64/N'=80 run used the N=48/N'=64 oee separator). It is required " *
                    "here because the Lehmann-Behnke path needs several dense interval copies and " *
                    "only ~159 GB is available on this host."),
        "quadrature"=>Dict("nt"=>48,"n_p"=>24,"npanel"=>npan,
            "npanel_rule"=>"npanel = round(4.8*LX); reproduces the published defaults 48 at LX=10 and 96 at LX=20",
            "t_star"=>ts,"t_star_rule"=>"t_star = 1.0 raised in 0.25 steps until the large-t moment tail bound <= 1e-13",
            "moment_tail_bound_at_t_star"=>tail,"clearance"=>clearance),
        "confinement"=>Dict("sigma"=>sig,"lambda_ref"=>c.lambda_ref,
            "margin"=>sig-c.lambda_ref,"holds"=>sig > c.lambda_ref,
            "formula"=>"sigma = inf_{x notin Omega} V, evaluated at the six face centres"),
        "reference"=>Dict("value"=>c.lambda_ref,"error"=>c.lambda_ref_err,"kind"=>c.ref_kind),
        "purpose"=>c.purpose,
        "goerisch_shift"=>Dict("c"=>CSHIFT(),
            "stageA_inf_L1"=>inf(L1),
            "certified_lambda_min_lower_bound"=>CSHIFT()+inf(L1),
            "rule"=>("c = max(1.0, ceil((-inf(L1)+0.25)*4)/4). H+cI must be CERTIFIED positive definite, "*
                     "i.e. c + inf(L1) > 0. c=1.0 is an H2+-tuned value that FAILS for H3^2+ at d_paper=2.8 "*
                     "(c+inf(L1) = -0.0515). c is a free parameter of the Goerisch construction, so raising it "*
                     "to satisfy the prerequisite is required by the method and is not a tuning of the result."),
            "raised_above_H2plus_default"=>CSHIFT() > 1.0),
        # Relative accuracy, so the claim is calibrated to what is being asserted.
        # For H3^2+ the stated tolerance is 5% relative; lambda_ref there is the
        # PROVISIONAL floating-Gaussian value (a corrected extended-precision value
        # is being computed separately), but a shift of order 1e-5 is irrelevant at
        # this tolerance.
        "relative_accuracy"=>Dict(
            "abs_lambda_ref"=>abs(c.lambda_ref),
            "certified_width"=>W,
            "certified_relative_width"=>ok ? W/abs(c.lambda_ref) : nothing,
            "certified_relative_width_percent"=>ok ? 100*W/abs(c.lambda_ref) : nothing,
            "target_relative"=>c.system == "H3plus" ? 0.05 : nothing,
            "target_absolute_width"=>c.system == "H3plus" ? 0.05*abs(c.lambda_ref) : nothing,
            "factor_inside_target"=>(ok && c.system == "H3plus") ? (0.05*abs(c.lambda_ref))/W : nothing,
            "meets_target"=>(ok && c.system == "H3plus") ? (W <= 0.05*abs(c.lambda_ref)) : nothing,
            "reference_status"=>c.system == "H3plus" ?
                "PROVISIONAL floating-Gaussian value; relative figures must be requoted if the corrected extended-precision value differs" :
                "exact"),
        "float_comparison"=>Dict(
            "source"=>"optimal_box_surface.csv, same (system, d_paper, L, N, N') cell",
            "float_LG_lower"=>c.float_lower, "float_dirichlet_upper"=>c.float_upper,
            "float_width"=>c.float_width,
            "certified_width"=>W,
            "interval_inflation_width_ratio"=>(ok && isfinite(c.float_width)) ? W/c.float_width : nothing,
            "certified_contains_float_lower"=>(ok && isfinite(c.float_lower)) ?
                (lam_lo <= c.float_lower <= lam_hi) : nothing,
            "certified_contains_float_upper"=>(ok && isfinite(c.float_upper)) ?
                (lam_lo <= c.float_upper <= lam_hi) : nothing,
            "lower_endpoint_gap_certified_minus_float"=>(ok && isfinite(c.float_lower)) ?
                lam_lo - c.float_lower : nothing,
            "upper_endpoint_gap_certified_minus_float"=>(ok && isfinite(c.float_upper)) ?
                lam_hi - c.float_upper : nothing),
        "decisive_test"=>(isfinite(c.ref_var_upper) || startswith(tag,"D")) ? Dict{String,Any}(
            "question"=>("does the CERTIFIED lower bound respect an independently verified variational "*
                         "UPPER bound on lambda_1? A rigorous lower bound must lie strictly below it."),
            "variational_upper_bound"=>isfinite(c.ref_var_upper) ? c.ref_var_upper : nothing,
            "best_estimate"=>c.lambda_ref, "best_estimate_error"=>c.lambda_ref_err,
            "certified_lower"=>lam_lo,
            "signed_gap_certified_lower_minus_variational_upper"=>
                (ok && isfinite(c.ref_var_upper)) ? lam_lo - c.ref_var_upper : nothing,
            "signed_gap_certified_lower_minus_best_estimate"=>ok ? lam_lo - c.lambda_ref : nothing,
            "VIOLATES_variational_upper_bound"=>
                (ok && isfinite(c.ref_var_upper)) ? (lam_lo > c.ref_var_upper) : nothing,
            "VIOLATES_best_estimate"=>ok ? (lam_lo > c.lambda_ref) : nothing,
            "float_lower_signed_gap_to_variational_upper"=>
                (isfinite(c.float_lower) && isfinite(c.ref_var_upper)) ? c.float_lower - c.ref_var_upper : nothing,
            "certification_moved_lower_bound_by"=>
                (ok && isfinite(c.float_lower)) ? lam_lo - c.float_lower : nothing,
            "note"=>("a POSITIVE signed gap means the certified lower bound EXCEEDS a verified variational "*
                     "upper bound, which is a correctness failure in the chain and not a width issue")) : nothing,
        "three_centre_note"=> c.system == "H3plus" ?
            "Three collinear unit charges at (-d,0,0), (0,0,0), (+d,0,0). No new mathematics " *
            "and no special branch was required: the centre at the ORIGIN is the case s=0 of the " *
            "SAME moments_verified.jl routine used for H2+ (moment_smallt integrates " *
            "cos(kappa(x+L))exp(-t^2(x-s)^2) over [-L,L] with its per-panel Bernstein remainder, " *
            "and moment_larget is the infinite-domain form plus a Gaussian tail bound; neither " *
            "has an s-dependent branch). The Coulomb singularity is never evaluated pointwise " *
            "because the Gaussian representation 1/r = (2/sqrt(pi)) int exp(-t^2 r^2) dt removes " *
            "it analytically. The only structural addition is the grouping of centres sharing a " *
            "transverse position, which sums their Bx factors before the Kronecker product and " *
            "keeps the cost identical to the two-centre case." : nothing,
        "CERTIFIED_ENCLOSURE"=>Dict(
            "lambda1_lower"=>lam_lo, "lambda1_upper"=>lam_hi, "width"=>W,
            "valid"=>ok,
            "lower_source"=>"$(prefix)_stagelg.json : lambda1_lower_certified (eee Lehmann-Goerisch, separator $(recL["separator_source"]))",
            "upper_source"=>"$(prefix)_stagedir.json : lambda1_upper_certified (Dirichlet Rayleigh quotient sup)",
            "contains_reference"=>contains_ref),
        "stages"=>Dict("A"=>recA,"mu2"=>recM,"lg"=>recL,"dirichlet"=>recD),
        "driver"=>Dict("script"=>"run_cert.jl","invocation"=>"julia -t <threads> run_cert.jl $tag",
            "core"=>"cert_core.jl","moments"=>"moments_verified.jl (bundle, unmodified)",
            "output_file"=>"$(prefix).json",
            "stage_files"=>Dict("A"=>"$(prefix)_stageA.json","mu2"=>"$(prefix)_stagemu2.json",
                               "lg"=>"$(prefix)_stagelg.json","dirichlet"=>"$(prefix)_stagedir.json")),
        "corrected_goerisch"=>Dict{String,Any}(
            "realization"=>"X = (L^2)^3 x L^2 x H^1(Omega); Tu = {sqrt(1-eps')grad u, sqrt(c-Ceps')u, u}; D = H^1(Omega)",
            "ceps_file"=>GCEPSFILE, "ceps_record"=>GCEPSREC,
            "eps_prime_grid"=>eps_grid_for_goerisch(),
            "discretization_term_included"=>GDISC,
            "gate_G1_mode"=>!GDISC,
            "uncorrected_A2_used"=>GOFF,
            "stageA_dirichlet_source"=>(GFROM == "" ? "recomputed in this run" : "read from $(GFROM)_stage{A,dir}.json"),
            "G4_self_consistency"=>G4),
        "threads"=>Threads.nthreads(), "julia"=>string(VERSION),
        "wall_seconds_total"=>time()-T0, "peak_rss_gb"=>rss())
    open("$(prefix).json","w") do f; JSON.print(f, jsan(out), 2); end
    @printf("\nwrote %s.json   total %.0f s  peak %.1f GB\n", prefix, out["wall_seconds_total"], rss())
    println("CERT_DONE_$(tag)")
    flush(stdout)
end

main()
