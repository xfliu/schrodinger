# moments_verified.jl — Track A / Stage-2, step A1
# Rigorous interval enclosure of the 1-D shifted cosine Coulomb moment
#     G(κ,t) = ∫_{-L}^{L} cos(κ (x+L)) e^{-t² (x-s)²} dx
# using ONLY exp/cos/sqrt (IntervalArithmetic 1.0.8) + two analytic error terms.
#
# Two regimes split at t_star:
#   large t (t ≥ t_star): infinite-domain closed form ± Gaussian tail bound
#   small t (t <  t_star): COMPOSITE verified Gauss–Legendre on [-L,L]
#                          (npanel panels, order n_p) + per-panel Bernstein remainder
#
# GL nodes/weights certified by interval-Newton on Legendre roots at BigFloat
# precision (endpoint-safe), then outward-rounded to Float64 intervals.
module MomentsVerified

using IntervalArithmetic: Interval, interval, mid, inf, sup, radius, diam,
                          intersect_interval, isempty_interval, issubset_interval

const IV  = Interval{Float64}
const IVB = Interval{BigFloat}

@inline PI() = interval(pi)
@inline SQRTPI() = sqrt(PI())

# ---- BigFloat point/interval Legendre (P_n and P_n') ----------------------
function _legp_big(n::Int, x::BigFloat)
    n==0 && return (BigFloat(1), BigFloat(0)); n==1 && return (x, BigFloat(1))
    p0=BigFloat(1); p1=x
    for k in 1:(n-1); p2=((2k+1)*x*p1 - k*p0)/(k+1); p0=p1; p1=p2; end
    dp = n*(x*p1 - p0)/(x*x - 1); (p1, dp)
end
function _legpd_ivb(n::Int, x::IVB)
    n==0 && return (interval(BigFloat(1)), interval(BigFloat(0)))
    n==1 && return (x, interval(BigFloat(1)))
    p0=interval(BigFloat(1)); p1=x
    for k in 1:(n-1)
        p2=(interval(BigFloat(2k+1))*x*p1 - interval(BigFloat(k))*p0)/interval(BigFloat(k+1))
        p0=p1; p1=p2
    end
    dp = interval(BigFloat(n))*(x*p1 - p0)/(x*x - interval(BigFloat(1))); (p1, dp)
end

# outward Float64 conversion of a BigFloat interval
_to_f64(Xb::IVB) = interval(Float64(inf(Xb), RoundDown), Float64(sup(Xb), RoundUp))

# ---- Certified reference GL rule on [-1,1] (order n) ----------------------
# returns Float64-interval nodes & weights, each provably enclosing the true value
function gl_reference(n::Int; prec::Int=512)
    setprecision(BigFloat, prec)
    nodes = Vector{IV}(undef, n); wts = Vector{IV}(undef, n)
    for i in 1:n
        x = BigFloat(cos(pi*(i-0.25)/(n+0.5)))
        for _ in 1:300
            p,dp = _legp_big(n, x); dx = p/dp; x -= dx
            abs(dx) < BigFloat(1e-70) && break
        end
        hw = BigFloat(1e-55); X = interval(x-hw, x+hw); ok=false
        for _ in 1:80
            m = interval(mid(X)); pm,_ = _legpd_ivb(n, m); _,dpX = _legpd_ivb(n, X)
            N = m - pm/dpX; Xn = intersect_interval(N, X)
            isempty_interval(Xn) && break
            if issubset_interval(N, X); X = Xn; ok=true; break; end
            X = diam(Xn) < diam(X) ? Xn : X
        end
        ok || error("gl_reference: interval-Newton failed to certify root $i (n=$n)")
        _,dp = _legpd_ivb(n, X)
        wb = interval(BigFloat(2))/((interval(BigFloat(1)) - X*X)*dp*dp)
        nodes[i] = _to_f64(X); wts[i] = _to_f64(wb)
    end
    return nodes, wts
end

# ---- Per-panel Bernstein remainder ----------------------------------------
# integrand cos(κ(x+L)) e^{-t²(x-s)²}, entire. On a panel of half-length h,
# ellipse E_ρ (in reference variable) gives η_max = h(ρ-1/ρ)/2.
# |f| ≤ e^{t² η²} · cosh(κ η) =: M(ρ). Trefethen: |E_n| ≤ 64 h M(ρ)/(15 ρ^{2n}(ρ²-1)).
icosh(a::IV) = (exp(a) + exp(-a)) / interval(2.0)
function panel_remainder(κ::Float64, t::Float64, h::Float64, n::Int)
    best=Inf; bestρ=1.5; ρ=1.05
    while ρ < 60.0
        b=(ρ-1/ρ)/2; η=h*b
        lb = log(64.0*h/15.0) + t*t*η*η + κ*η - 2n*log(ρ) - log(ρ*ρ-1)
        if lb < best; best=lb; bestρ=ρ; end
        ρ *= 1.02
    end
    ρI = interval(bestρ); bI=(ρI - interval(1.0)/ρI)/interval(2.0); ηI=interval(h)*bI
    MI = exp(interval(t*t)*ηI*ηI) * icosh(interval(κ)*ηI)
    EI = interval(64.0)*interval(h)*MI/(interval(15.0)*ρI^(2n)*(ρI*ρI-interval(1.0)))
    E = sup(EI); return interval(-E, E)
end

# ---- Small-t: composite verified GL + Bernstein remainder -----------------
function moment_smallt(κ::Float64, t::Float64, L::Float64, s::Float64,
                       nodes::Vector{IV}, wts::Vector{IV}, npanel::Int)
    n_p = length(nodes)
    h  = L/npanel                          # panel half-length
    hI = interval(h); κI=interval(κ); tI=interval(t); sI=interval(s)
    acc = interval(0.0)
    for j in 0:(npanel-1)
        cj = -L + (2j+1)*h                 # panel centre (float; enclosed below)
        cjI = interval(cj)
        for i in 1:n_p
            xI = cjI + hI*nodes[i]
            fI = cos(κI*(xI + interval(L))) * exp(-tI*tI*(xI - sI)*(xI - sI))
            acc += (hI * wts[i]) * fI
        end
    end
    acc += interval(npanel) * panel_remainder(κ, t, h, n_p)   # total remainder ≤ npanel·max
    return acc
end

# ---- Large-t: infinite-domain closed form ± Gaussian tail -----------------
function moment_larget(κ::Float64, t::Float64, L::Float64, s::Float64)
    LI=interval(L); κI=interval(κ); tI=interval(t); sI=interval(s)
    Ginf = cos(κI*(sI+LI)) * (SQRTPI()/tI) * exp(-(κI*κI)/(interval(4.0)*tI*tI))
    dL = LI+sI; dR = LI-sI
    tail = (interval(1.0)/(tI*tI*dL))*exp(-tI*tI*dL*dL) +
           (interval(1.0)/(tI*tI*dR))*exp(-tI*tI*dR*dR)
    tmag = sup(tail); return Ginf + interval(-tmag, tmag)
end

# ---- Dispatcher -----------------------------------------------------------
function moment_verified(κ::Float64, t::Float64, L::Float64, s::Float64;
                         t_star::Float64=1.0,
                         nodes::Union{Nothing,Vector{IV}}=nothing,
                         wts::Union{Nothing,Vector{IV}}=nothing,
                         npanel::Int=96)
    if t >= t_star
        return moment_larget(κ, t, L, s)
    else
        (nodes===nothing || wts===nothing) && error("moment_verified: small-t needs certified nodes/wts")
        return moment_smallt(κ, t, L, s, nodes, wts, npanel)
    end
end

end # module
