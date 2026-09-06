# run_formfactor.jl -- CERTIFIED tail-refined form factor eta' (item 5 / P1).
#
#   julia -t <threads> run_formfactor.jl <N> <sigma> <sector:eee|oee> <LX> <nt> [mu list]
#
# DEFINITION USED, AND THE ASSUMPTION IN IT.  The paper's Remark app:rem:si-tail is not in the
# reproduction bundle, so Pi^0_N could not be checked against it.  This file assumes
# Pi^0_N = the L^2-ORTHOGONAL projection onto V_N, which is what the L^2 norm in the numerator
# ||(I - Pi^0_N) V v_N|| implies.  With an ORTHONORMAL cosine basis (M = I) that gives
#     ||(I - Pi^0_N) V v_N||^2 = ||V v_N||^2 - ||Pi_N(V v_N)||^2 = x'(P2 - P^2)x
# so, with a_sigma = a_0 + sigma(.,.) having Gram matrix S = K + P + sigma I,
#     eta'^2 = lambda_max of the pencil (Q', S),   Q' = P2 - P^2.
# If instead Pi^0_N is the a_0-orthogonal projection, the middle factor is not I and the
# (2,2) block below changes from I to the a_0 Gram matrix -- the SAME cost, one extra
# Cholesky, so it can be added on request.  FLAGGED, not silently assumed away.
#
# CERTIFICATION ROUTE (a rigorous UPPER bound, which is the conservative direction).
# lambda_max(Q',S) <= mu  iff  mu S - Q' = mu S - P2 + P^2 >= 0.  Forming P^2 is impossible at
# this size (D^3 interval FMA), so use the Schur complement: since I > 0,
#     [[ mu S - P2 ,  P ] ,  [ P , I ]]  >= 0   <==>   mu S - P2 - P I^{-1} P >= 0,
# i.e. exactly the condition wanted -- and the 2D x 2D block matrix NEVER forms P^2.
# PSD of the resulting symmetric INTERVAL matrix [A] = Am +- Ar is verified by
#     lambda_min(Am) > ||Ar||_2,     ||Ar||_2 <= max_i sum_j Ar_ij   (Ar >= 0 entrywise),
# and lambda_min(Am) > rho is certified by a FLOAT Cholesky of Am - (rho + c) I completing with
# all pivots positive, where c = (n+1) u max_i |Am_ii| / (1 - (n+1) u), u = 2^-53, is the
# standard backward-error allowance for Cholesky (Rump, verification of positive definiteness).
# Only float Cholesky is needed, so the cost is one O((2D)^3/3) factorisation per trial mu.
using LinearAlgebra, Printf, JSON
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam, radius
push!(LOAD_PATH, get(ENV, "VEIGS_SRC", ""))
include("cert_core.jl"); using .CertCore; const CC = CertCore
const IV = Interval{Float64}
rss() = CC.peak_rss_gb()

# ---------------------------------------------------------- the V^2 MATRIX ----
# Same collapse as the streamed form, but accumulated into a D x D interval matrix because
# the pencil needs Q' as an operator, not a scalar quadratic form.
function assemble_P2(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                     centres, LX, LY, LZ, nt::Int, npanel::Int, t_star::Float64,
                     nodes, wts, verbose::Bool=true)
    @assert LY == LZ "fast path assumes LY == LZ"
    t, wt = CC.build_t_grid(nt)
    Ix, Iy, Iz = CC.sector_idx(N,px), CC.sector_idx(N,py), CC.sector_idx(N,pz)
    nx, ny, nz = length(Ix), length(Iy), length(Iz)
    D = nx*ny*nz
    pref2 = (interval(2.0)/sqrt(interval(pi)))^2
    ax = [c[1][1] for c in centres]; Zs = [c[2] for c in centres]
    P2 = fill(interval(0.0), D, D)
    npair = 0; nterm = 0; npairs_tot = nt*(nt+1)÷2
    for k in 1:nt, l in k:nt
        t2 = t[k]^2; u2 = t[l]^2; w2 = t2 + u2; wid = sqrt(w2); gam = t2*u2/w2
        base = interval((k == l) ? 1.0 : 2.0) * pref2 *
               interval(wt[k]) * interval(wt[l])
        By = CC.cos_moment_matrix_1d(N, LY, 0.0, wid, nodes, wts;
                                     npanel=npanel, t_star=t_star)[Iy,Iy]
        Bz = By          # LY == LZ asserted, and py == pz, so Iy == Iz
        # In ANY single-parity axis both |n-m| and n+m are EVEN, so only even f enters G[f]
        # and the x -> -x substitution costs sin(f pi/2) = 0: B(centre s) = B(centre -s).
        buckets = Dict{Float64,IV}()
        for i in eachindex(ax), j in eachindex(ax)
            m = (t2*ax[i] + u2*ax[j])/w2
            wij = interval(Zs[i]*Zs[j]) * exp(-interval(gam)*interval((ax[i]-ax[j])^2))
            key = abs(m); buckets[key] = get(buckets, key, interval(0.0)) + wij
        end
        for (m, wij) in buckets
            Bx = CC.cos_moment_matrix_1d(N, LX, m, wid, nodes, wts;
                                         npanel=npanel, t_star=t_star)[Ix,Ix]
            CC.kron_accum_iv!(P2, base*wij, Bx, By, Bz)
            nterm += 1
        end
        npair += 1
        verbose && npair % 100 == 0 &&
            (@printf("      [P2] pair %d/%d terms %d peak %.1f GB\n",
                     npair, npairs_tot, nterm, rss()); flush(stdout))
    end
    return P2, D, nterm, npair
end

# ---------------------------------------------- verified PSD of the Schur block ---
# Builds Am = mid([[mu S - P2, P],[P, I]]) and the radius row-sums, then Cholesky.
# Pm, P2m, Kv are float midpoints; Pr, P2r are float radii.
# CORRECTED ROUTE.  The condition is mu S - Q' = mu S - P2 + P^2 >= 0.  A Schur complement
# SUBTRACTS its off-diagonal square, so the block [[mu S - P2, P],[P, I]] encodes
# mu S - P2 - P^2 -- the WRONG SIGN, and it is unverifiable for any mu (measured: NO_MU_VERIFIED
# even at mu = 1 where the true value is 0.35).  Instead form P^2 as a FLOAT gemm and carry a
# rigorous perturbation bound, which is also cheaper: dimension D, not 2D.
#
# For any P in [P], P2 in [P2], with E1 = P - Pm, E2 = P2 - P2m and Chat = fl(Pm*Pm):
#   A_true = mu(K + P + sigma I) - P2 + P^2
#   A_true - Ahat = mu E1 - E2 + (Pm E1 + E1 Pm + E1^2) + (Pm^2 - Chat)
# so, using ||X||_2 <= max_i sum_j |X_ij| for each nonnegative bound matrix,
#   ||A_true - Ahat||_2 <= mu r1 + r2 + 2 pn r1 + r1^2 + D u pn^2 =: rho_tot,
# with r1 = max row sum of rad(P), r2 = max row sum of rad(P2), pn = max row sum of |Pm|,
# and D u pn^2 the standard gemm rounding bound.  Then all members are PSD once
# lambda_min(Ahat) > rho_tot, certified by float Cholesky of Ahat - (rho_tot + c)I completing.
function verify_mu(mu::Float64, Pm::Matrix{Float64}, P2m::Matrix{Float64}, Kv::Vector{Float64},
                   Pr::Vector{Float64}, P2r::Vector{Float64}, sigma::Float64,
                   pn::Float64, Chat::Matrix{Float64})
    D = size(Pm,1)
    A = Matrix{Float64}(undef, D, D)
    Threads.@threads for j in 1:D
        @inbounds for i in 1:D
            A[i,j] = mu*Pm[i,j] - P2m[i,j] + Chat[i,j]
        end
        @inbounds A[j,j] += mu*(Kv[j] + sigma)
    end
    r1 = maximum(Pr); r2 = maximum(P2r)
    u = 2.0^-53
    rho = mu*r1 + r2 + 2*pn*r1 + r1*r1 + D*u*pn*pn
    dmax = 0.0
    @inbounds for i in 1:D; dmax = max(dmax, abs(A[i,i])); end
    nu = (D+1)*u
    cc = nu*dmax/(1.0 - nu)
    shift = rho + cc
    @inbounds for i in 1:D; A[i,i] -= shift; end
    ok = true
    try
        cholesky!(Symmetric(A, :L))
    catch e
        ok = false
    end
    A = nothing; GC.gc(true)
    return ok, rho, cc, shift, dmax
end

function main()
    T0 = time()
    N     = parse(Int, ARGS[1]); sigma = parse(Float64, ARGS[2])
    sect  = ARGS[3];             LX    = parse(Float64, ARGS[4])
    nt    = parse(Int, ARGS[5])
    mus   = length(ARGS) >= 6 ? [parse(Float64,x) for x in split(ARGS[6],",")] :
            [0.2080, 0.2100, 0.2150, 0.2250, 0.2500, 0.3000]
    LY = 0.8*LX; LZ = LY
    px, py, pz = sect == "eee" ? (:even,:even,:even) : (:odd,:even,:even)
    cen  = CC.h2plus_centres(2.0)
    npan = Int(round(4.8*LX)); ts, tail = CC.choose_t_star(cen, LX, LY, LZ)
    nodes, wts = CC.MomentsVerified.gl_reference(24)
    @printf("===== certified tail-refined eta' : N=%d sigma=%.4f sector=%s LX=%.1f nt=%d =====\n",
            N, sigma, sect, LX, nt)
    @printf("  npanel=%d t_star=%.2f cell_tail=%.2e threads=%d BLAS=%d julia=%s\n",
            npan, ts, tail, Threads.nthreads(), BLAS.get_num_threads(), VERSION); flush(stdout)

    t1 = @elapsed (P, K, D) = CC.assemble_PK(N,px,py,pz; centres=cen, LX=LX, LY=LY, LZ=LZ,
                                 npanel=npan, t_star=ts, nodes=nodes, wts=wts)
    @printf("  [V matrix] D=%d (%.1f s, peak %.1f GB)\n", D, t1, rss()); flush(stdout)
    # certified Galerkin mu_1, mu_2 of this sector: eigen of the midpoint plus interval radius
    Hm = Matrix{Float64}(undef, D, D)
    Threads.@threads for j in 1:D
        @inbounds for i in 1:D; Hm[i,j] = mid(P[i,j]); end
        @inbounds Hm[j,j] += mid(K[j])
    end
    Hr = Vector{Float64}(undef, D)
    Threads.@threads for i in 1:D
        s = 0.0; @inbounds for j in 1:D; s += radius(P[i,j]); end
        Hr[i] = s
    end
    radH = maximum(Hr)
    ev = eigen(Symmetric(Hm), 1:2)
    mu1 = ev.values[1] - radH; mu2 = ev.values[2] - radH   # rigorous lower bounds
    mu1u = ev.values[1] + radH
    @printf("  [Galerkin] mu_1^N in [%.14f, %.14f], mu_2^N >= %.14f (||rad||_2 <= %.3e)\n",
            mu1, mu1u, mu2, radH); flush(stdout)
    Hm = nothing; GC.gc(true)

    # float midpoints/radii of P, then free the interval matrix
    Pm = Matrix{Float64}(undef, D, D); Pr = Vector{Float64}(undef, D)
    Threads.@threads for j in 1:D
        @inbounds for i in 1:D; Pm[i,j] = mid(P[i,j]); end
    end
    Threads.@threads for i in 1:D
        s = 0.0; @inbounds for j in 1:D; s += radius(P[i,j]); end
        Pr[i] = s
    end
    Kv = [mid(K[j]) for j in 1:D]
    P = nothing; K = nothing; GC.gc(true)
    @printf("  [P -> float] peak %.1f GB\n", rss()); flush(stdout)

    t2 = @elapsed (P2, D2, nterm, npair) = assemble_P2(N,px,py,pz; centres=cen, LX=LX, LY=LY,
            LZ=LZ, nt=nt, npanel=npan, t_star=ts, nodes=nodes, wts=wts)
    @assert D2 == D
    @printf("  [V^2 matrix] %d terms over %d node pairs (%.1f s, peak %.1f GB)\n",
            nterm, npair, t2, rss()); flush(stdout)
    P2m = Matrix{Float64}(undef, D, D); P2r = Vector{Float64}(undef, D)
    Threads.@threads for j in 1:D
        @inbounds for i in 1:D; P2m[i,j] = mid(P2[i,j]); end
    end
    Threads.@threads for i in 1:D
        s = 0.0; @inbounds for j in 1:D; s += radius(P2[i,j]); end
        P2r[i] = s
    end
    P2 = nothing; GC.gc(true)
    @printf("  [P2 -> float] peak %.1f GB\n", rss()); flush(stdout)

    # P^2 as a float gemm (BLAS-3), plus the row-sum norm of |Pm| for the perturbation bound
    pn = 0.0
    let rows = Vector{Float64}(undef, D)
        Threads.@threads for i in 1:D
            s_ = 0.0; @inbounds for j in 1:D; s_ += abs(Pm[i,j]); end
            rows[i] = s_
        end
        pn = maximum(rows)
    end
    tg = @elapsed Chat = Pm*Pm
    @printf("  [P^2 float gemm] %.1f s, ||Pm||_rowsum = %.6f, peak %.1f GB\n", tg, pn, rss())
    flush(stdout)

    # ---- direct float lambda_max(Q',S) whenever forming P^2 densely is affordable ----
    direct = nothing
    if D <= 9000
        Qf = P2m - Chat
        Sf = copy(Pm); @inbounds for j in 1:D; Sf[j,j] += Kv[j] + sigma; end
        dv = eigen(Symmetric(Qf), Symmetric(Sf)).values
        direct = maximum(dv)
        @printf("  [direct float] lambda_max(Q',S) = %.12f  -> eta' = %.10f  (D=%d, dense P^2)\n",
                direct, sqrt(max(direct,0.0)), D); flush(stdout)
        Qf = nothing; Sf = nothing; GC.gc(true)
    end

    # ------------------------------------------------- the certified upper bound ---
    trials = Dict{String,Any}(); best = nothing
    for mu in mus
        tv = @elapsed ((ok, rho, cc, shift, dmax) = verify_mu(mu, Pm, P2m, Kv, Pr, P2r, sigma, pn, Chat))
        @printf("  mu=%.6f : PSD verified = %-5s  (radius bound %.3e, allowance %.3e, %.1f s, peak %.1f GB)\n",
                mu, ok, rho, cc, tv, rss()); flush(stdout)
        trials["mu=$mu"] = Dict("verified"=>ok, "radius_bound_2norm"=>rho,
            "cholesky_allowance"=>cc, "total_diag_shift"=>shift, "max_abs_diag"=>dmax,
            "seconds"=>tv)
        if ok; best = mu; break; end
    end

    # ---- feed eta'^2 through Theorem 8.  THE THEOREM'S EXACT STATEMENT IS NOT IN THE BUNDLE.
    # The existing code path (CertCore.sharp_L_iv) builds its form factor as
    #   etaV^2 = eps/(1-eps) + Ceps/ms1,  g = etaV/sqrt(ms1),
    #   Ch2 = (1+g^2)/((1-eps) nu_*(N)),  L_k = msk/(1+Ch2 msk) - shift.
    # Substituting the CERTIFIED eta'^2 for etaV^2 leaves one ambiguity: whether the refined
    # constant keeps the (1-eps) denominator, which in the analytic route came from the eps
    # splitting that the direct measurement no longer uses.  BOTH variants are reported.
    thm8 = Dict{String,Any}()
    if best !== nothing
        nus = ((interval(N+1)*interval(pi))/(interval(2.0)*interval(LX)))^2
        for (vn, fac) in (("no_eps_factor", interval(1.0)), ("keeps_1_minus_eps_0p4", interval(0.6)))
            ms1 = interval(mu1) + interval(sigma); ms2 = interval(mu2) + interval(sigma)
            g2  = interval(best)/ms1
            Ch2 = (interval(1.0) + g2)/(fac*nus)
            L1  = ms1/(interval(1.0) + Ch2*ms1) - interval(sigma)
            L2  = ms2/(interval(1.0) + Ch2*ms2) - interval(sigma)
            thm8[vn] = Dict("Ch2"=>sup(Ch2), "L1_lower"=>inf(L1), "L2_lower"=>inf(L2))
            @printf("  [Theorem 8 / %s] Ch2=%.8e  L_1 >= %.12f  L_2 >= %.12f\n",
                    vn, sup(Ch2), inf(L1), inf(L2)); flush(stdout)
        end
        thm8["CAVEAT"] = "the paper's Theorem 8 statement is NOT in the reproduction bundle. These "*
            "two rows are the EXISTING code path (CertCore.sharp_L_iv) with the certified eta'^2 "*
            "substituted for its analytic etaV^2, under the two readings of the (1-eps) factor. "*
            "Neither is claimed to be the paper's Theorem 8; compare against the independent "*
            "verifier's L_1 = -0.673784 / L_2 = -0.476988 before use."
        thm8["verifier_reference"] = Dict("eta_prime"=>0.4555, "L1"=>-0.673784, "L2"=>-0.476988,
                                          "U1"=>-0.551851)
    end

    out = Dict{String,Any}(
      "theorem8_evaluation"=>thm8,
      "N"=>N, "sigma"=>sigma, "sector"=>sect, "LX"=>LX, "box"=>[LX,LY,LZ], "D"=>D, "nt"=>nt,
      "quadrature"=>Dict("npanel"=>npan,"t_star"=>ts,"cell_tail"=>tail,"n_p"=>24),
      "definition"=>Dict(
        "eta_prime"=>"sup_{v_N in V_N} ||(I - Pi^0_N) V v_N|| / ||v_N||_{a_sigma}",
        "Pi_0_N_ASSUMED"=>"L^2-orthogonal projection onto V_N (the paper's Remark app:rem:si-tail "*
                          "is NOT in the reproduction bundle, so this could not be checked)",
        "consequence"=>"eta'^2 = lambda_max(Q', S), Q' = P2 - P^2, S = K + P + sigma I, M = I "*
                       "because the cosine basis is orthonormal",
        "if_a0_projection"=>"the (2,2) block of the Schur matrix changes from I to the a_0 Gram "*
                            "matrix; same cost, one extra Cholesky, not done here"),
      "certification_route"=>Dict(
        "statement"=>"lambda_max(Q',S) <= mu  iff  mu S - P2 + P^2 >= 0. P^2 is formed as a FLOAT "*
                     "gemm and the whole perturbation is bounded rigorously, so no interval "*
                     "matrix product is needed and the dimension stays D.",
        "perturbation_bound"=>"||A_true - Ahat||_2 <= mu r1 + r2 + 2 pn r1 + r1^2 + D u pn^2, with "*
                     "r1 = max row sum of rad(P), r2 = max row sum of rad(P2), pn = max row sum of "*
                     "|Pm|, from A_true - Ahat = mu E1 - E2 + (Pm E1 + E1 Pm + E1^2) + (Pm^2 - Chat)",
        "psd_test"=>"lambda_min(Ahat) > rho_tot certified by float Cholesky of "*
                    "Ahat - (rho_tot + c)I completing with positive pivots, "*
                    "c = (D+1)u max|Ahat_ii|/(1-(D+1)u), u = 2^-53",
        "rejected_route"=>"the Schur block [[mu S - P2, P],[P, I]] encodes mu S - P2 - P^2 -- the "*
                     "Schur complement SUBTRACTS, so the sign of P^2 is wrong. Measured: it is "*
                     "unverifiable for every mu up to 1.0 at N=16 and N=24, where the true "*
                     "lambda_max is 0.35. Recorded so the mistake is not repeated.",
        "direction"=>"an UPPER bound on eta'^2, i.e. conservative in the safe direction"),
      "galerkin"=>Dict("mu1_lower"=>mu1, "mu1_upper"=>mu1u, "mu2_lower"=>mu2,
                       "interval_radius_2norm_bound"=>radH),
      "trials"=>trials,
      "direct_float_lambda_max"=>direct,
      "direct_float_eta_prime"=>direct === nothing ? nothing : sqrt(max(direct,0.0)),
      "direct_note"=>"forms Q' = P2 - P^2 densely in FLOAT and solves the generalised eigenproblem; "*
                     "only affordable for D <= 9000, and is the cross-check that the certified "*
                     "upper bound brackets the true value from above",
      "eta_prime_sq_certified_upper"=>best,
      "eta_prime_certified_upper"=>best === nothing ? nothing : sqrt(best),
      "timings"=>Dict("V_matrix"=>t1,"V2_matrix"=>t2,"total"=>time()-T0),
      "peak_rss_gb"=>rss(),
      "host"=>Dict("julia"=>string(VERSION),"threads"=>Threads.nthreads()))
    fn = "formfactor_$(sect)_N$(N)_sig$(replace(string(sigma),"."=>"p")).json"
    open(fn,"w") do f; JSON.print(f, out, 2); end
    @printf("\nwrote %s  total %.0f s  peak %.1f GB\n", fn, out["timings"]["total"], rss())
    best === nothing ? println("NO_MU_VERIFIED") :
        @printf("CERTIFIED eta'^2 <= %.10f  (eta' <= %.10f)\n", best, sqrt(best))
    println("FORMFACTOR_DONE_$(sect)_$(N)")
end
main()
