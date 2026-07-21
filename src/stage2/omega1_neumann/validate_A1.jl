
using Printf, LinearAlgebra
using SpecialFunctions: erfcx
include("moments_verified.jl")
using .MomentsVerified
using IntervalArithmetic: inf, sup, diam, interval

function Gfloat(κ,t,L,s)
    α=-L-s; β=L-s
    za=t*α-im*κ/(2t); basea=exp(-t^2*α^2+im*κ*α)
    Ba= α>=0 ? erfcx(za)*basea : 2*exp(-κ^2/(4t^2))-erfcx(-za)*basea
    zb=t*β-im*κ/(2t); baseb=exp(-t^2*β^2+im*κ*β)
    Bb= β>=0 ? erfcx(zb)*baseb : 2*exp(-κ^2/(4t^2))-erfcx(-zb)*baseb
    real(exp(im*κ*(L+s))*(sqrt(pi)/(2t))*(Ba-Bb))
end
function build_t_grid(nt)
    k=1:(nt-1); β=collect(k ./ sqrt.(4.0 .* k.^2 .- 1.0))
    J=SymTridiagonal(zeros(nt),β); v,V=eigen(J); x=v; w=2.0 .* (V[1,:].^2)
    u=0.5 .*(x .+1); wu=0.5 .*w; (u ./(1 .-u), wu ./(1 .-u).^2)
end

L=10.0; t_all,_=build_t_grid(48); t_star=1.0
n_p=24; npanel=48
tc=@elapsed nodes,wts=MomentsVerified.gl_reference(n_p)
@printf("GL ref n_p=%d: maxnodew=%.2e maxwtw=%.2e (%.2fs)\n", n_p, maximum(diam.(nodes)), maximum(diam.(wts)), tc)

# full grid containment at N=32 and N=64 (frequencies up to 2N)
for N in (32, 64)
    fmax=2N; nfail=0; ntest=0; maxw=0.0; worst=(0.0,0.0); maxsmall=0.0; maxlarge=0.0
    tt=@elapsed for t in t_all, f in 0:fmax, s in (-2.0,2.0)
        κ=f*pi/(2L)
        G=MomentsVerified.moment_verified(κ,t,L,s; t_star=t_star, nodes=nodes, wts=wts, npanel=npanel)
        g=Gfloat(κ,t,L,s); ntest+=1
        (inf(G)<=g<=sup(G)) || (nfail+=1)
        w=diam(G); if w>maxw; maxw=w; worst=(κ,t) end
        if t>=t_star; maxlarge=max(maxlarge,w) else maxsmall=max(maxsmall,w) end
    end
    @printf("N=%d fmax=%d tests=%d FAIL=%d maxwidth=%.3e (κ=%.2f t=%.2e) small=%.3e large=%.3e (%.1fs)\n",
            N,fmax,ntest,nfail,maxw,worst[1],worst[2],maxsmall,maxlarge,tt)
end
println("A1_DONE")
