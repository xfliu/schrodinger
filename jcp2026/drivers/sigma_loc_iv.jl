# sigma_loc_iv.jl -- rigorous closed-form coercivity constant, in interval arithmetic.
#
# REPLACES the previous shift, which took C_eps from the auxiliary Rayleigh-Ritz value
#     Ceps[e] = interval(-sup(eta), -inf(eta)),  eta = lambda_1(e*K + P)
# That quantity is INADMISSIBLE: the cosine space is a subspace of H^1(Omega), so
# Rayleigh-Ritz bounds from ABOVE, eta_{1,N} >= eta, hence -eta_{1,N} <= C_eps^opt -- it sits
# strictly BELOW the constant the coercivity argument requires, at every N.  The auxiliary
# eigenvalue is retained downstream as a DIAGNOSTIC ONLY, under the key
# "eta_1N_ritz_diagnostic_not_a_bound".
#
# Lemma (manuscript).  For nuclei a_i with dist(a_i, boundary) >= rho, any 0 < r1 < r2 <= rho
# and any delta > 0,
#
#     (V^- v, v) <= eps' |grad v|^2 + sigma_loc |v|^2      for all v in H^1(Omega),
#
#     sigma_loc(eps') = Ztot^2 (1+delta)/(4 eps')
#                     + (1 + 1/delta) eps' / ((1+delta) (r2-r1)^2)
#                     + Ztot / r1 .
#
# Valid for ANY admissible (delta, r1, r2) -- no optimisation is needed for validity, and
# none is performed here.  The parameters are fixed by the caller and the closed form is
# evaluated in interval arithmetic.
module SigmaLoc

using IntervalArithmetic: Interval, interval, inf, sup, diam
const IV = Interval{Float64}

export sigma_loc_iv, box_rho, total_charge

"""
    box_rho(centres, LX, LY, LZ) -> Float64

rho = min over nuclei of the distance from that nucleus to the boundary of the box
[-LX,LX] x [-LY,LY] x [-LZ,LZ].  Computed from the geometry record; never hardcoded.
For Omega_1 ([-10,10]x[-8,8]^2, nuclei at (+-2,0,0)) this returns min(8,8,8) = 8;
for Omega_2 ([-20,20]x[-16,16]^2) it returns min(18,16,16) = 16.
"""
function box_rho(centres, LX::Float64, LY::Float64, LZ::Float64)
    r = Inf
    for (a, _) in centres
        r = min(r, LX - abs(a[1]), LY - abs(a[2]), LZ - abs(a[3]))
    end
    return r
end

total_charge(centres) = sum(z for (_, z) in centres)

"""
    sigma_loc_iv(eps, delta, r1, r2, rho, Ztot) -> IV

Rigorous enclosure of sigma_loc.  The three hypotheses of the lemma are asserted; the caller
is expected to let the run die rather than proceed with an inadmissible parameter set.

`delta` may be passed as a Rational to get an exact interval for the intended rational value
(e.g. 1//10) rather than the nearest Float64.
"""
function sigma_loc_iv(eps::Real, delta::Real, r1::Real, r2::Real, rho::Real, Ztot::Real)
    # --- hypotheses of the lemma, asserted at run time -----------------------
    @assert eps   > 0            "sigma_loc hypothesis violated: eps must be > 0 (got $eps)"
    @assert delta > 0            "sigma_loc hypothesis violated: delta must be > 0 (got $delta)"
    @assert r1    > 0            "sigma_loc hypothesis violated: need 0 < r1 (got r1=$r1)"
    @assert r1    < r2           "sigma_loc hypothesis violated: need r1 < r2 (got r1=$r1, r2=$r2)"
    @assert r2    <= rho         "sigma_loc hypothesis violated: need r2 <= rho (got r2=$r2, rho=$rho)"
    @assert Ztot  > 0            "sigma_loc hypothesis violated: Ztot must be > 0 (got $Ztot)"

    iv(x::Rational) = interval(numerator(x)) / interval(denominator(x))
    iv(x::Real)     = interval(Float64(x))

    e  = iv(eps);  d  = iv(delta)
    R1 = iv(r1);   R2 = iv(r2);   Z = iv(Ztot)
    one_ = interval(1.0);  four = interval(4.0)

    gap = R2 - R1
    t1 = Z*Z*(one_ + d)/(four*e)                       # Ztot^2 (1+delta)/(4 eps)
    t2 = (one_ + one_/d)*e/((one_ + d)*gap*gap)        # (1+1/delta) eps / ((1+delta)(r2-r1)^2)
    t3 = Z/R1                                          # Ztot / r1
    s  = t1 + t2 + t3

    # Free self-check: (1+1/d)/(1+d) == 1/d identically, so t2 must also equal eps/(d*gap^2).
    # If the two interval evaluations were disjoint something is badly wrong with the
    # arithmetic; they are only required to overlap, since each carries its own rounding.
    t2alt = e/(d*gap*gap)
    @assert !(sup(t2) < inf(t2alt) || sup(t2alt) < inf(t2)) "sigma_loc internal check failed: " *
        "t2=[$(inf(t2)),$(sup(t2))] disjoint from equivalent form [$(inf(t2alt)),$(sup(t2alt))]"

    return s
end

"""
    sigma_loc_terms(eps, delta, r1, r2, rho, Ztot) -> Dict

Per-term breakdown for the output record, so the constant is auditable rather than a
single opaque number.
"""
function sigma_loc_terms(eps::Real, delta::Real, r1::Real, r2::Real, rho::Real, Ztot::Real)
    iv(x::Rational) = interval(numerator(x)) / interval(denominator(x))
    iv(x::Real)     = interval(Float64(x))
    e = iv(eps); d = iv(delta); R1 = iv(r1); R2 = iv(r2); Z = iv(Ztot)
    one_ = interval(1.0); four = interval(4.0); gap = R2 - R1
    t1 = Z*Z*(one_ + d)/(four*e)
    t2 = (one_ + one_/d)*e/((one_ + d)*gap*gap)
    t3 = Z/R1
    s  = sigma_loc_iv(eps, delta, r1, r2, rho, Ztot)
    return Dict{String,Any}(
        "sigma_loc"       => [inf(s), sup(s)],
        "sigma_loc_width" => diam(s),
        "term_gradient_Ztot2_1pd_over_4eps" => [inf(t1), sup(t1)],
        "term_cutoff_1p1od_eps_over_1pd_gap2" => [inf(t2), sup(t2)],
        "term_far_field_Ztot_over_r1" => [inf(t3), sup(t3)],
        "parameters" => Dict{String,Any}("eps"=>Float64(eps), "delta"=>Float64(delta),
            "delta_exact"=>string(delta), "r1"=>Float64(r1), "r2"=>Float64(r2),
            "rho"=>Float64(rho), "Ztot"=>Float64(Ztot)),
        "hypotheses" => Dict{String,Any}("eps_gt_0"=>eps>0, "delta_gt_0"=>delta>0,
            "0_lt_r1"=>r1>0, "r1_lt_r2"=>r1<r2, "r2_le_rho"=>r2<=rho),
        "status" => "RIGOROUS closed-form coercivity constant (manuscript lemma); " *
                    "valid for any admissible (delta,r1,r2), no optimisation performed")
end

end # module
