# run_sweep.jl -- (L,N) enclosure-width sweep driver.
#
#   julia -t <threads> run_sweep.jl <spec.csv> <out_prefix>
#
# spec.csv lines:  system,d_paper,L,N,Nprime,tier
#   system : H2plus | H3plus
#   tier   : dense  -> Neumann eee mu1 Ritz (+ eee 2nd) and Dirichlet eee Ritz
#            lg     -> dense + oee mu2 Ritz + full LG lower at aux N'
#            lg_eoe -> lg + eoe sector ground (separator-sector spot check)
#
# Writes <out_prefix>_rows.csv (one row per cell, appended and flushed as each
# cell finishes, so a crash costs only the cell in flight) and
# <out_prefix>_detail.jsonl (per-cell detail: rho sweep, A0/A2, iteration counts).
using LinearAlgebra, Printf
include("sweep_core.jl")
using .SweepCore
const SC = SweepCore
const FC = SweepCore.FloatCore

# ---------------------------------------------------------------- systems ----
# Reference lambda_1(R^3) in PAPER units, per system/geometry.  PROVENANCE:
#   H2plus  : exact prolate-spheroidal 1ssg at standard R=2 bohr, -1.1026342144949,
#             halved by the paper convention -> -0.5513170.  Source:
#             h3plus_reference.json .validation_H2plus.reference_1ssg_std, itself
#             from the project artifact h2plus_R2_bounds.csv.
#   H3plus  : floating s-type GTO Rayleigh-Ritz values from the artifact
#             h3plus_reference.json .study_geometries[].lambda1_paper with
#             .error_estimate_paper.  There is NO published benchmark for this
#             system, so these ARE the reference; each is a variational UPPER
#             bound, so the true value lies in [value - err, value], and `err`
#             is the floor on any agreement claim made against it.
function lambda_ref(system::String, d::Float64)
    if system == "H2plus"
        return -0.5513170, 0.0            # exact spheroidal, half of -1.1026342
    elseif system == "H3plus"
        if isapprox(d, 4.0; atol=1e-9);  return -0.76207999734,  1.5e-7; end
        if isapprox(d, 2.8; atol=1e-9);  return -0.936459672195, 5.0e-8; end
        error("no reference lambda_1 for H3plus d_paper=$d")
    end
    error("unknown system $system")
end

centres_of(system::String, d::Float64) =
    system == "H2plus" ? SC.h2plus_centres(2.0) : SC.h3plus_centres(d)

# separator allowance: the deficit actually realised by the published certified
# mu2 separator on Omega2 at N=64 --
#   mu2N_Ritz(oee) = -0.3333398828  vs  certified mu2 lower = -0.3395656252
const DELTA_SEP = 0.0062257424
const CSHIFT    = 1.0

# ------------------------------------------------------------- tiny JSON -----
jq(x::Float64) = isfinite(x) ? @sprintf("%.17g", x) : "null"
jq(x::Int) = string(x); jq(x::Bool) = x ? "true" : "false"
jq(x::String) = "\"" * replace(x, "\"" => "'") * "\""
jq(x::Vector) = "[" * join(jq.(x), ",") * "]"
jq(d::Dict) = "{" * join(["\"$k\":" * jq(d[k]) for k in sort(collect(keys(d)))], ",") * "}"

# ------------------------------------------------------- primary block -------
function primary_block_indices(N::Int, Nprime::Int)
    nxs = length(FC.sector_idx(N, :even)); nxb = length(FC.sector_idx(Nprime, :even))
    idx = Vector{Int}(undef, nxs^3); c = 0
    for cz in 1:nxs, cy in 1:nxs, cx in 1:nxs
        c += 1; idx[c] = ((cz-1)*nxb + (cy-1))*nxb + cx
    end
    return idx, nxs, nxb
end

const HDR = ["system","d_paper","box","L","Ly","N","Nprime","D","D_aux","npanel","t_star",
             "mu1N_Ritz","eee2_Ritz","lambda1D_Ritz","mu2N_oee_Ritz","eoe_Ritz",
             "rho_used","LG_lower","LG_lower_ideal_rho","enclosure_width",
             "truncation_term","LG_discretisation_term","trunc_gap_discrete",
             "lambda_ref","sigma_confinement","hypothesis_holds","max_a_over_L",
             "moment_tail_bound","tier","status","reason","wall_seconds"]

blank(x) = x === nothing ? "" :
           x isa Float64 ? @sprintf("%.17g", x) :
           x isa Bool    ? (x ? "true" : "false") :
           x isa String  ? (occursin(',', x) ? "\"" * x * "\"" : x) : string(x)

function main()
    spec_path = ARGS[1]; prefix = ARGS[2]
    nb = haskey(ENV,"SWEEP_BLAS") ? parse(Int, ENV["SWEEP_BLAS"]) : 48
    NT = haskey(ENV,"SWEEP_NT")   ? parse(Int, ENV["SWEEP_NT"])   : 48
    @printf("outer t-grid nt=%d\n", NT)
    BLAS.set_num_threads(nb)
    @printf("BLAS threads=%d\n", BLAS.get_num_threads())
    nodes, wts = FC.gl_rule(24)
    rows_path = prefix * "_rows.csv"; det_path = prefix * "_detail.jsonl"
    open(rows_path, "w") do f; println(f, join(HDR, ",")); end
    open(det_path, "w") do f; end

    cells = Tuple{String,Float64,Float64,Int,Int,String}[]
    for ln in eachline(spec_path)
        (isempty(strip(ln)) || startswith(ln, "#") || startswith(ln, "system")) && continue
        p = split(strip(ln), ",")
        push!(cells, (String(p[1]), parse(Float64,p[2]), parse(Float64,p[3]),
                      parse(Int,p[4]), parse(Int,p[5]), String(p[6])))
    end
    @printf("threads=%d julia=%s cells=%d\n", Threads.nthreads(), VERSION, length(cells))
    flush(stdout)

    for (ci,(system,d,L,N,Np,tier)) in enumerate(cells)
        T0 = time()
        LX = L; LY = 0.8*L; LZ = 0.8*L
        cen = centres_of(system, d)
        lref, lerr = lambda_ref(system, d)
        npan = SC.npanel_rule(LX)
        ts, tail = SC.choose_t_star(cen, LX, LY, LZ; tol=1e-13)
        maxaL = maximum(abs(c[1][1]) for c in cen) / LX
        sig = sigma_conf(cen, LX, LY, LZ)
        holds = sig > lref
        row = Dict{String,Any}("system"=>system,"d_paper"=>d,
            "box"=>@sprintf("[-%g,%g]x[-%g,%g]^2", LX,LX,LY,LY),
            "L"=>LX,"Ly"=>LY,"N"=>N,"Nprime"=>nothing,"D"=>nothing,"D_aux"=>nothing,
            "npanel"=>npan,"t_star"=>ts,"lambda_ref"=>lref,
            "sigma_confinement"=>sig,"hypothesis_holds"=>holds,
            "max_a_over_L"=>maxaL,"moment_tail_bound"=>tail,"tier"=>tier,
            "status"=>"ok","reason"=>"")
        for k in ("mu1N_Ritz","eee2_Ritz","lambda1D_Ritz","mu2N_oee_Ritz","eoe_Ritz",
                  "rho_used","LG_lower","LG_lower_ideal_rho","enclosure_width",
                  "truncation_term","LG_discretisation_term","trunc_gap_discrete")
            row[k] = nothing
        end
        det = Dict{String,Any}("system"=>system,"d_paper"=>d,"L"=>LX,"N"=>N,"tier"=>tier)

        @printf("\n[%d/%d] %s d=%.3g L=%g N=%d tier=%s npanel=%d t*=%.2f tail=%.2e sigma=%.6f\n",
                ci, length(cells), system, d, LX, N, tier, npan, ts, tail, sig); flush(stdout)

        if !holds
            row["status"] = "skipped"
            row["reason"] = @sprintf("confinement fails: sigma=%.6f <= lambda_ref=%.6f", sig, lref)
            @printf("   SKIPPED: %s\n", row["reason"]); flush(stdout)
        else
          try
            shift0 = lref - 0.35
            # ---------------- Neumann eee ----------------
            ta = @elapsed (H, D) = SC.assemble_neumann_mc(N,:even,:even,:even; centres=cen,
                    LX=LX,LY=LY,LZ=LZ,npanel=npan,t_star=ts,nodes=nodes,wts=wts,nt=NT)
            tv = @elapsed (vals,_,it,inf1) = SC.ground_si(H, shift0; nev=2)
            H = nothing; GC.gc()
            mu1 = vals[1]; eee2 = vals[2]
            row["D"] = D; row["mu1N_Ritz"] = mu1; row["eee2_Ritz"] = eee2
            det["neumann"] = Dict{String,Any}("D"=>D,"t_assemble"=>ta,"t_eig"=>tv,
                "iters"=>it,"converged"=>inf1.converged,"shift"=>shift0)
            @printf("   mu1N=%.13f eee2=%.6f  (asm %.1fs eig %.1fs it=%d conv=%s)\n",
                    mu1, eee2, ta, tv, it, inf1.converged); flush(stdout)

            # ---------------- Dirichlet eee ----------------
            tad = @elapsed (HD, DD, nxd) = SC.assemble_dirichlet_mc(N; centres=cen,
                    LX=LX,LY=LY,LZ=LZ,npanel=npan,t_star=ts,nodes=nodes,wts=wts,nt=NT)
            tvd = @elapsed (dvals,_,itd,infd) = SC.ground_si(HD, shift0; nev=1)
            HD = nothing; GC.gc()
            l1d = dvals[1]
            row["lambda1D_Ritz"] = l1d
            row["trunc_gap_discrete"] = l1d - mu1
            det["dirichlet"] = Dict{String,Any}("D"=>DD,"modes_per_axis"=>nxd,
                "t_assemble"=>tad,"t_eig"=>tvd,"iters"=>itd,"converged"=>infd.converged)
            @printf("   l1D =%.13f  gap=%+.4e  (asm %.1fs eig %.1fs it=%d)\n",
                    l1d, l1d-mu1, tad, tvd, itd, ); flush(stdout)

            if tier == "lg" || tier == "lg_eoe"
                # ---------------- oee mu2 (separator source) ----------------
                tao = @elapsed (HO, DO) = SC.assemble_neumann_mc(N,:odd,:even,:even; centres=cen,
                        LX=LX,LY=LY,LZ=LZ,npanel=npan,t_star=ts,nodes=nodes,wts=wts,nt=NT)
                tvo = @elapsed (ovals,_,ito,_) = SC.ground_si(HO, shift0; nev=1)
                HO = nothing; GC.gc()
                mu2 = ovals[1]; row["mu2N_oee_Ritz"] = mu2
                det["oee"] = Dict{String,Any}("D"=>DO,"t_assemble"=>tao,"t_eig"=>tvo,"iters"=>ito)
                @printf("   mu2N(oee)=%.13f  (asm %.1fs eig %.1fs)\n", mu2, tao, tvo); flush(stdout)

                if tier == "lg_eoe"
                    tae = @elapsed (HE, DE) = SC.assemble_neumann_mc(N,:even,:odd,:even; centres=cen,
                            LX=LX,LY=LY,LZ=LZ,npanel=npan,t_star=ts,nodes=nodes,wts=wts,nt=NT)
                    (evals,_,ite,_) = SC.ground_si(HE, shift0; nev=1)
                    HE = nothing; GC.gc()
                    row["eoe_Ritz"] = evals[1]
                    det["eoe"] = Dict{String,Any}("D"=>DE,"t_assemble"=>tae,"ground"=>evals[1])
                    @printf("   eoe ground=%.13f  (oee<eoe: %s)\n", evals[1], mu2 < evals[1])
                    flush(stdout)
                end

                # ---------------- LG lower at aux N' ----------------
                taa = @elapsed (Hc, Db) = SC.assemble_neumann_mc(Np,:even,:even,:even; centres=cen,
                        LX=LX,LY=LY,LZ=LZ,npanel=npan,t_star=ts,nodes=nodes,wts=wts,nt=NT,
                        sigma=CSHIFT)
                row["Nprime"] = Np; row["D_aux"] = Db
                @printf("   aux N'=%d D=%d assembled (%.1fs)\n", Np, Db, taa); flush(stdout)
                pidx, nxs, nxb = primary_block_indices(N, Np)
                Hsub = Hc[pidx,pidx]
                tps = @elapsed (pvals,pvecs,itp,_) = SC.ground_si(Hsub, shift0+CSHIFT; nev=1)
                Hsub = nothing; GC.gc()
                v = zeros(Float64, Db); v[pidx] .= pvecs[:,1]; v ./= norm(v)
                A0 = dot(v, Symmetric(Hc)*v)
                tcg = @elapsed (w, nit) = FC.cg_solve(Hc, v; reltol=1e-12, maxiter=8000)
                cgres = norm(Symmetric(Hc)*w .- v)/norm(v)
                Hw = Symmetric(Hc)*w
                A2 = 2.0*dot(w,v) - dot(w,Hw)
                Hc = nothing; GC.gc()
                rho_ideal = mu2
                rho_cert  = mu2 - DELTA_SEP
                Li,_,Bi,nui = FC.lg_mobius(A0, A2, rho_ideal, CSHIFT)
                Lc,_,Bc,nuc = FC.lg_mobius(A0, A2, rho_cert,  CSHIFT)
                row["rho_used"] = rho_cert
                row["LG_lower"] = Lc; row["LG_lower_ideal_rho"] = Li
                row["enclosure_width"] = l1d - Lc
                row["truncation_term"] = l1d - lref
                row["LG_discretisation_term"] = lref - Lc
                rsweep = Dict{String,Any}()
                for f in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0)
                    r = mu2 - f*DELTA_SEP
                    Lr,_,Br,nur = FC.lg_mobius(A0, A2, r, CSHIFT)
                    rsweep[@sprintf("%.2f", f)] = [r, Lr, Br, nur]
                end
                det["lg"] = Dict{String,Any}("Nprime"=>Np,"D_aux"=>Db,"A0"=>A0,"A2"=>A2,
                    "cg_iters"=>nit,"cg_residual"=>cgres,"t_aux_assemble"=>taa,
                    "t_primary_eig"=>tps,"t_cg"=>tcg,"primary_iters"=>itp,
                    "rho_ideal"=>rho_ideal,"rho_cert"=>rho_cert,
                    "LG_ideal"=>Li,"LG_cert"=>Lc,"B_ideal"=>Bi,"B_cert"=>Bc,
                    "nu_ideal"=>nui,"nu_cert"=>nuc,
                    "B_positive_cert"=>(Bc>0),"nu_lt_1_cert"=>(nuc<1),
                    "rho_sweep_factor_of_DELTA_SEP"=>rsweep)
                @printf("   LG: A0=%.12f A2=%.12f cg=%d/%.1e | rho_ideal=%.10f -> %.13f | rho_cert=%.10f -> %.13f\n",
                        A0, A2, nit, cgres, rho_ideal, Li, rho_cert, Lc)
                @printf("   W = l1D - LG_cert = %.6e   (trunc %.3e + LGdisc %.3e)\n",
                        l1d-Lc, l1d-lref, lref-Lc); flush(stdout)
            end
          catch e
            row["status"] = "failed"
            row["reason"] = first(replace(string(typeof(e)) * ": " * sprint(showerror, e), "," => ";", "\n" => " "), 220)
            @printf("   FAILED: %s\n", row["reason"]); flush(stdout)
          end
        end
        row["wall_seconds"] = time() - T0
        det["wall_seconds"] = row["wall_seconds"]; det["status"] = row["status"]
        open(rows_path, "a") do f
            println(f, join([blank(row[h]) for h in HDR], ","))
        end
        open(det_path, "a") do f; println(f, jq(det)); end
    end
    println("SWEEP_DONE")
    flush(stdout)
end

# sigma(Omega) = inf_{x not in Omega} V, attained on the boundary.  V is smooth
# on each open face and the nuclei are all on the x-axis, so the minimum over a
# transverse (y or z) face is at its centre and the minimum over an x-face is at
# its centre too; sigma is the smaller of the two.  Validated against the
# recorded H2+ values -0.242535625 (Omega1) and -0.124034735 (Omega2).
function sigma_conf(cen, LX, LY, LZ)
    V(x,y,z) = -sum(c[2]/sqrt((x-c[1][1])^2 + (y-c[1][2])^2 + (z-c[1][3])^2) for c in cen)
    return min(V(LX,0.0,0.0), V(-LX,0.0,0.0), V(0.0,LY,0.0), V(0.0,-LY,0.0),
               V(0.0,0.0,LZ), V(0.0,0.0,-LZ))
end

main()
