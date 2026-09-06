# run_formfactor_h3_sect.jl -- CERTIFIED tail-refined form factor eta' (item 5 / P1).
#
#   julia -t <threads> run_formfactor_h3_sect.jl <N> <sigma> <sector:eee|oee|eoe|eeo> <LX> <nt> [mu list]
#
# COPY of run_formfactor_h3.jl (md5 5ab31d38a0a84d0bbfea8f2c1694f1cd) with FOUR changes, all
# marked in place, made so that the MIXED transverse sectors eoe / eeo can be run:
#   (i)   sector parsing accepts eee | oee | eoe | eeo and ERRORS on anything else (the
#         original silently ran oee for every string other than "eee");
#   (ii)  the `Bz = By` copy in BOTH V^2 assemblies (interval assemble_P2 and float
#         assemble_P2_float) is now taken only under an explicit `py === pz` guard.  LY == LZ
#         was @assert-ed but py == pz was not, so a mixed sector silently used the y-parity
#         block on the z axis.  Under the guard the shortcut returns the SAME object, so the
#         eee/oee arithmetic -- and its bits -- are unchanged (gate G1);
#   (iii) FF_GEOM gained `h3z:<dp>` (three centres, ALL charges zero) for the V = 0 exactness
#         gate G2, together with a t_star guard for the all-zero case (CertCore.cell_tail
#         reduces over the nonzero charges and throws on an empty collection);
#   (iv)  FF_NEV (default 0) adds a SEPARATE eigen(Symmetric(Hm), 1:FF_NEV) and records the
#         enclosures under "galerkin_extra_ritz".  With FF_NEV unset nothing is computed and
#         no key is added, so the certificate has exactly the original's field set.
# NOT changed: the bucketing on abs(m).  It merges the +s and -s Gaussian centres and its
# justification needs the axis carrying the centres to be single-parity; the centres differ
# only along x and every sector above is EVEN in x, so the original argument holds verbatim.
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
using Dates
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
        # (ii) THE Bz FIX.  The original wrote `Bz = By` under the comment "LY == LZ asserted,
        # and py == pz, so Iy == Iz".  LY == LZ is @assert-ed above but py == pz was NEVER
        # checked, so for a mixed transverse sector (eoe: py=:odd, pz=:even) the y-parity block
        # was silently used on the z axis.  The copy shortcut is now taken ONLY under an
        # explicit py === pz guard, so eee/oee cost and BITS are unchanged.
        Bz = py === pz ? By :
             CC.cos_moment_matrix_1d(N, LZ, 0.0, wid, nodes, wts;
                                     npanel=npanel, t_star=t_star)[Iz,Iz]
        # In ANY single-parity axis both |n-m| and n+m are EVEN, so only even f enters G[f]
        # and the x -> -x substitution costs sin(f pi/2) = 0: B(centre s) = B(centre -s).
        buckets = Dict{Float64,IV}()
        for i in eachindex(ax), j in eachindex(ax)
            Zs[i]*Zs[j] == 0.0 && continue     # a zero charge contributes nothing (bundle rule)
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


# ------------------------------------- FLOAT Q assembly with a-priori rounding ---
# Q is assembled in Float64 from the SAME interval-certified 1-D moment matrices, and the
# floating-point error is covered by the standard a-priori bound.  For n accumulated terms
#   |fl(Q) - Q|_ij <= gamma_n * (sum_s |c_s| |Bx_s| |By| |Bz|)_ij + (radius contributions)
# with gamma_n = n u / (1 - n u).  Both bounds are accumulated as SCALARS via elementwise
# maxima -- conservative, and it keeps memory at ONE D x D Float64 matrix instead of three.
# Then ||E||_2 <= ||E||_F <= D * max_ij |E_ij|.
# The s-quadrature enclosure is untouched; only the D x D accumulation changes arithmetic.
function assemble_P2_float(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                           centres, LX, LY, LZ, nt::Int, npanel::Int, t_star::Float64,
                           nodes, wts, verbose::Bool=true)
    @assert LY == LZ "fast path assumes LY == LZ"
    t, wt = CC.build_t_grid(nt)
    Ix, Iy, Iz = CC.sector_idx(N,px), CC.sector_idx(N,py), CC.sector_idx(N,pz)
    nx, ny, nz = length(Ix), length(Iy), length(Iz)
    D = nx*ny*nz
    # index maps, matching the column-major (x fastest) ordering of kron_accum_iv!
    JX = Vector{Int}(undef, D); JY = Vector{Int}(undef, D); JZ = Vector{Int}(undef, D)
    @inbounds for i in 1:D
        q, r = divrem(i-1, nx); z, y = divrem(q, ny)
        JX[i] = r+1; JY[i] = y+1; JZ[i] = z+1
    end
    pref2 = (interval(2.0)/sqrt(interval(pi)))^2
    ax = [c[1][1] for c in centres]; Zs = [c[2] for c in centres]
    Q = zeros(Float64, D, D)
    sabs = 0.0; srad = 0.0
    npair = 0; nterm = 0; npairs_tot = nt*(nt+1)÷2
    for k in 1:nt, l in k:nt
        t2 = t[k]^2; u2 = t[l]^2; w2 = t2 + u2; wid = sqrt(w2); gam = t2*u2/w2
        base = interval((k == l) ? 1.0 : 2.0) * pref2 *
               interval(wt[k]) * interval(wt[l])
        Byi = CC.cos_moment_matrix_1d(N, LY, 0.0, wid, nodes, wts;
                                      npanel=npanel, t_star=t_star)[Iy,Iy]
        Bym = mid.(Byi); byr = maximum(radius.(Byi)); bym = maximum(abs.(Bym))
        # (ii) THE Bz FIX, float path.  The original reused Bym on the z axis as well; that is
        # only legitimate when py == pz (and LY == LZ, @assert-ed above).  Under the guard the
        # copy is the SAME object, so every float operation below is bitwise unchanged for
        # eee/oee; for a mixed sector the z block is built from its own parity index set.
        Bzi = py === pz ? Byi :
              CC.cos_moment_matrix_1d(N, LZ, 0.0, wid, nodes, wts;
                                      npanel=npanel, t_star=t_star)[Iz,Iz]
        Bzm = py === pz ? Bym : mid.(Bzi)
        bzr = py === pz ? byr : maximum(radius.(Bzi))
        bzm = py === pz ? bym : maximum(abs.(Bzm))
        buckets = Dict{Float64,IV}()
        for i in eachindex(ax), j in eachindex(ax)
            Zs[i]*Zs[j] == 0.0 && continue     # a zero charge contributes nothing (bundle rule)
            m = (t2*ax[i] + u2*ax[j])/w2
            wij = interval(Zs[i]*Zs[j]) * exp(-interval(gam)*interval((ax[i]-ax[j])^2))
            key = abs(m); buckets[key] = get(buckets, key, interval(0.0)) + wij
        end
        for (m, wij) in buckets
            Bxi = CC.cos_moment_matrix_1d(N, LX, m, wid, nodes, wts;
                                          npanel=npanel, t_star=t_star)[Ix,Ix]
            ci = base*wij
            c = mid(ci); cr = radius(ci)
            Bxm = mid.(Bxi); bxr = maximum(radius.(Bxi)); bxm = maximum(abs.(Bxm))
            Threads.@threads for j in 1:D
                jx = JX[j]; jy = JY[j]; jz = JZ[j]
                @inbounds for i in 1:D
                    Q[i,j] += c * Bxm[JX[i],jx] * Bym[JY[i],jy] * Bzm[JZ[i],jz]
                end
            end
            sabs += abs(c)*bxm*bym*bzm
            srad += cr*bxm*bym*bzm +
                    abs(c)*(bxr*bym*bzm + bxm*byr*bzm + bxm*bym*bzr)
            nterm += 1
        end
        npair += 1
        verbose && npair % 100 == 0 &&
            (@printf("      [Qf] pair %d/%d terms %d peak %.1f GB\n",
                     npair, npairs_tot, nterm, rss()); flush(stdout))
    end
    uu = 2.0^-53
    gam_n = nterm*uu/(1.0 - nterm*uu)
    Emax = gam_n*sabs + srad
    @printf("  [Qf bound] nterm=%d gamma_n=%.4e sabs=%.4e srad=%.4e -> max|E|=%.4e, ||E||_2<=%.4e\n",
            nterm, gam_n, sabs, srad, Emax, D*Emax)
    flush(stdout)
    return Q, D, nterm, npair, Emax
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
    # (i) SECTOR PARSING.  The original admitted only "eee" and let EVERY other string fall
    # through to oee.  All four sectors with at most one odd direction are now named, and an
    # unrecognised string is an error instead of a silent oee run.
    px, py, pz = sect == "eee" ? (:even,:even,:even) :
                 sect == "oee" ? (:odd, :even,:even) :
                 sect == "eoe" ? (:even,:odd, :even) :
                 sect == "eeo" ? (:even,:even,:odd ) :
                 error("unrecognised sector \"$sect\": expected one of eee, oee, eoe, eeo")
    # ---- geometry selection (FF_GEOM).  The default reproduces the promoted driver EXACTLY.
    #   h2:<ax>    two unit charges at x = -ax, +ax                       (default "h2:2.0")
    #   h3:<dp>    linear H3^2+: three unit charges at x = -dp, 0, +dp     (paper units)
    #   h3z0:<ax>  REGRESSION geometry: h2:<ax> plus a ZERO charge at the origin
    geom = get(ENV, "FF_GEOM", "h2:2.0")
    gkind, gpar = split(geom, ":"); gval = parse(Float64, gpar)
    cen = gkind == "h2"   ? CC.h2plus_centres(gval) :
          gkind == "h3"   ? CC.h3plus_centres(gval) :
          gkind == "h3z0" ? CC.Centre[((-gval,0.0,0.0),1.0), ((0.0,0.0,0.0),0.0), ((gval,0.0,0.0),1.0)] :
          gkind == "h3z"  ? CC.Centre[((-gval,0.0,0.0),0.0), ((0.0,0.0,0.0),0.0), ((gval,0.0,0.0),0.0)] :
          error("FF_GEOM must be h2:<ax>, h3:<dp>, h3z0:<ax> or h3z:<dp>; got $geom")
    gtag = geom == "h2:2.0" ? "" : "_" * replace(replace(geom, ":" => "d"), "." => "p")
    @printf("  [geometry] FF_GEOM=%s centres=%s tag='%s'\n", geom, string(cen), gtag); flush(stdout)
    # G2 (V = 0) support: CertCore.cell_tail reduces over `c for c in cen if c[2] != 0.0`, which
    # throws on an empty collection when EVERY charge is zero.  There is no potential to bound in
    # that case, so t_star keeps its starting value and the tail is exactly 0.  Any geometry with
    # a nonzero charge takes the original branch, so this is a no-op for eee/oee/production runs.
    allzero = all(cc -> cc[2] == 0.0, cen)
    npan = Int(round(4.8*LX))
    ts, tail = allzero ? (1.0, 0.0) : CC.choose_t_star(cen, LX, LY, LZ)
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
    # OPTIONAL extra Ritz enclosures (G2 needs three).  A SEPARATE eigen call placed AFTER the
    # original one, so the code path that produces mu1/mu1u/mu2 is untouched; and the JSON key
    # below is added ONLY when FF_NEV is set, so with FF_NEV unset the certificate has exactly
    # the original's fields.
    nev_extra = parse(Int, get(ENV, "FF_NEV", "0"))
    ritz_extra = nothing
    if nev_extra > 0
        eex = eigen(Symmetric(Hm), 1:nev_extra)
        ritz_extra = [[eex.values[i] - radH, eex.values[i] + radH] for i in 1:nev_extra]
        for i in 1:nev_extra
            @printf("  [Ritz %d] enclosure [%.16e, %.16e]\n", i, ritz_extra[i][1], ritz_extra[i][2])
        end
        flush(stdout)
    end
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

    # ---- V^2 cache.  Q = (V^2 phi_i, phi_j) depends on (LX, N, sector, nt, npanel, t_star)
    # and NOT on sigma or the mu ladder, so one build serves every sigma at that geometry.
    cdir  = get(ENV, "FF_CACHE", @__DIR__)
    ckey  = "v2_LX$(LX)_N$(N)_$(sect)_nt$(nt)_np$(npan)_ts$(ts)$(gtag).bin"
    cfile = joinpath(cdir, ckey)
    usec  = get(ENV, "FF_NOCACHE", "0") != "1"
    local P2m, P2r, t2, nterm, npair
    if usec && isfile(cfile)
        t2 = @elapsed begin
            open(cfile, "r") do io
                Dc = read(io, Int); Dc == D || error("V^2 cache dimension $Dc != $D in $cfile")
                nterm = read(io, Int); npair = read(io, Int)
                P2m = Array{Float64}(undef, D, D); read!(io, P2m)
                P2r = Vector{Float64}(undef, D);   read!(io, P2r)
            end
        end
        @printf("  [V^2 CACHE HIT] %s (%.1f s, peak %.1f GB)\n", ckey, t2, rss())
        flush(stdout)
    else
        usef = get(ENV, "FF_FLOATQ", "0") == "1"
        if get(ENV, "FF_QCHECK", "0") == "1"
            # GATE: the float Q with its E-ball must CONTAIN the interval Q elementwise.
            (Qi, _, _, _)          = assemble_P2(N,px,py,pz; centres=cen, LX=LX, LY=LY, LZ=LZ,
                                        nt=nt, npanel=npan, t_star=ts, nodes=nodes, wts=wts,
                                        verbose=false)
            (Qc, _, nt2, _, Em)    = assemble_P2_float(N,px,py,pz; centres=cen, LX=LX, LY=LY,
                                        LZ=LZ, nt=nt, npanel=npan, t_star=ts, nodes=nodes,
                                        wts=wts, verbose=false)
            worst = 0.0; nbad = 0
            @inbounds for j in 1:D, i in 1:D
                dev = abs(Qc[i,j] - mid(Qi[i,j]))
                tol = Em + radius(Qi[i,j])
                r = dev/tol
                r > worst && (worst = r)
                r > 1.0 && (nbad += 1)
            end
            @printf("  [Qf GATE] max ratio = %.6e over %d entries, violations %d -> %s\n",
                    worst, D*D, nbad, nbad == 0 ? "PASS" : "FAIL")
            flush(stdout)
            Qi = nothing; Qc = nothing; GC.gc(true)
            nbad == 0 || error("float Q containment gate FAILED at $nbad entries")
        end
        if usef
            t2 = @elapsed (Qf, D2, nterm, npair, Emax) = assemble_P2_float(N,px,py,pz;
                    centres=cen, LX=LX, LY=LY, LZ=LZ, nt=nt, npanel=npan, t_star=ts,
                    nodes=nodes, wts=wts)
            @assert D2 == D
            P2m = Qf
            P2r = fill(D*Emax, D)          # row-sum bound: D entries each <= Emax
            @printf("  [Qf FLOAT] %d terms over %d node pairs (%.1f s, peak %.1f GB), row bound %.4e\n",
                    nterm, npair, t2, rss(), D*Emax)
            flush(stdout)
        else
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
        end
        if usec
            tmp = cfile * ".tmp"
            open(tmp, "w") do io
                write(io, D); write(io, nterm); write(io, npair)
                write(io, P2m); write(io, P2r)
            end
            mv(tmp, cfile; force=true)
            @printf("  [V^2 CACHE WRITE] %s (%.2f GB)\n", ckey, filesize(cfile)/2^30)
            flush(stdout)
        end
    end

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
      "geometry"=>Dict("FF_GEOM"=>geom, "tag"=>gtag,
                       "centres"=>[Dict("x"=>c[1][1],"y"=>c[1][2],"z"=>c[1][3],"Z"=>c[2]) for c in cen]),
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
    if ritz_extra !== nothing
        out["galerkin_extra_ritz"] = Dict("nev"=>nev_extra, "enclosures"=>ritz_extra,
            "note"=>"[val - radH, val + radH] for the first FF_NEV Ritz values of this block, "*
                    "from a SEPARATE eigen(Symmetric(Hm), 1:FF_NEV) call that does not touch "*
                    "the mu1/mu2 path above")
    end
    fn = "formfactor$(gtag)_$(sect)_N$(N)_sig$(replace(string(sigma),"."=>"p")).json"
    open(fn,"w") do f; JSON.print(f, out, 2); end
    # PER-RUN COPY: the live name is keyed only by (sector, N, sigma), so a repeat run at the
    # same configuration overwrites it -- and a cache-warm re-run would replace a cold run's
    # timings with a cache-hit figure.  Write an immutable timestamped copy alongside.
    stamp = Dates.format(Dates.now(Dates.UTC), "yyyymmddTHHMMSSZ")
    fnprc = replace(fn, r"\.json$" => "_" * stamp * ".json")
    open(fnprc,"w") do f; JSON.print(f, out, 2); end
    @printf("\nwrote %s  total %.0f s  peak %.1f GB\n", fn, out["timings"]["total"], rss())
    @printf("per-run copy: %s\n", fnprc)
    best === nothing ? println("NO_MU_VERIFIED") :
        @printf("CERTIFIED eta'^2 <= %.10f  (eta' <= %.10f)\n", best, sqrt(best))
    println("FORMFACTOR_DONE_$(sect)_$(N)")
end
main()
