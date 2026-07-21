# dirichlet_verified.jl — rigorous Dirichlet (sine-basis) Rayleigh–Ritz UPPER bound
# on the H2+ ground energy. Provides the upper side of the domain-truncation chain
#     [Neumann LG lower] <= mu1 <= lambda1(R^3) <= lambda1^D <= lambda1^D_N  [this].
#
# Dirichlet sine basis psi_n(x) = (1/sqrt(L)) sin(n*pi*(x+L)/(2L)), n=1,2,...  (vanish at x=+-L).
#   kinetic eigenvalue (n*pi/(2L))^2 ; L2-orthonormal (mass = identity).
#   Coulomb sine moment: 0.5*(1/L)*(G[|n-m|] - G[n+m])   (sin*sin product-to-sum: MINUS).
#   Even-in-x functions are the ODD sine modes n=1,3,5,... (parity of sin(n*pi/2+theta)).
# Rigorous upper bound = sup( <v,H v>_interval / <v,v>_interval ) for the float eee ground
#   vector v; v in H^1_0 => its Rayleigh quotient is a rigorous upper bound on lambda1^D >= lambda1.
module DirichletVerified
using LinearAlgebra, Printf
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
include("moments_verified.jl")
using .MomentsVerified
const IV = Interval{Float64}
const LX, LY, LZ = 10.0, 8.0, 8.0
const AX = 2.0

build_t_grid(nt::Int) = begin
    k=1:(nt-1); β=collect(k ./ sqrt.(4.0 .*k.^2 .-1.0))
    J=SymTridiagonal(zeros(nt),β); vals,vecs=eigen(J)
    x=vals; w=2.0 .*(vecs[1,:].^2); u=0.5 .*(x .+1); wu=0.5 .*w
    (u ./(1 .-u), wu ./(1 .-u).^2)
end

# sine moment matrix over modes 1..mmax at fixed t (interval)
function sine_moment_matrix_1d(mmax::Int, L::Float64, s::Float64, t::Float64,
                               nodes::Vector{IV}, wts::Vector{IV}; npanel::Int=48, t_star::Float64=1.0)
    fmax = 2*mmax
    G = Vector{IV}(undef, fmax+1)
    for f in 0:fmax
        κ = f*pi/(2L)
        G[f+1] = MomentsVerified.moment_verified(κ, t, L, s; t_star=t_star, nodes=nodes, wts=wts, npanel=npanel)
    end
    Ninv = interval(1.0)/interval(L)     # (1/sqrt(L))^2
    half = interval(0.5)
    B = Matrix{IV}(undef, mmax, mmax)
    @inbounds for n in 1:mmax, m in 1:mmax
        B[n,m] = half*Ninv*(G[abs(n-m)+1] - G[n+m+1])
    end
    return B
end

# odd sine modes (even-in-x) among 1..mmax
odd_modes(mmax::Int) = collect(1:2:mmax)

# assemble Dirichlet eee interval Hamiltonian (odd sine modes per axis)
function assemble_dirichlet(N::Int; nt::Int=48, npanel::Int=48, nodes::Vector{IV}, wts::Vector{IV})
    mmax = N+1                              # gives (N/2+1) odd modes => matches Neumann eee D
    t, wt = build_t_grid(nt)
    Ix = odd_modes(mmax); Iy = Ix; Iz = Ix
    nx=length(Ix); D=nx^3
    pref = interval(2.0)/sqrt(interval(pi))
    P = fill(interval(0.0), D, D)
    for k in 1:nt
        Bxm = sine_moment_matrix_1d(mmax, LX, -AX, t[k], nodes, wts; npanel=npanel, t_star=t_star_())
        Bxp = sine_moment_matrix_1d(mmax, LX, +AX, t[k], nodes, wts; npanel=npanel, t_star=t_star_())
        Bx = (Bxm .+ Bxp)[Ix,Ix]
        By = sine_moment_matrix_1d(mmax, LY, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star_())[Iy,Iy]
        Bz = sine_moment_matrix_1d(mmax, LZ, 0.0, t[k], nodes, wts; npanel=npanel, t_star=t_star_())[Iz,Iz]
        wk = -pref*interval(wt[k])
        P .+= wk .* kron(Bz, kron(By, Bx))
    end
    # kinetic diagonal (n*pi/(2L))^2 for the selected odd modes
    νx=(interval(pi)/(interval(2.0)*interval(LX)))^2
    νy=(interval(pi)/(interval(2.0)*interval(LY)))^2
    νz=(interval(pi)/(interval(2.0)*interval(LZ)))^2
    mx=Ix; my=Iy; mz=Iz
    H = copy(P)
    @inbounds for cz in 1:nx, cy in 1:nx, cx in 1:nx
        idx=(cz-1)*nx*nx+(cy-1)*nx+cx
        H[idx,idx] = H[idx,idx] + νx*interval(mx[cx]^2)+νy*interval(my[cy]^2)+νz*interval(mz[cz]^2)
    end
    return H, D, nx
end
t_star_() = 1.0

imatvec(M::Matrix{IV}, x::Vector{IV}) = begin
    D=length(x); y=Vector{IV}(undef,D)
    @inbounds for i in 1:D
        acc=interval(0.0); @simd for j in 1:D; acc+=M[i,j]*x[j]; end; y[i]=acc
    end; y
end
idot(a::Vector{IV},b::Vector{IV})=(s=interval(0.0);@inbounds for i in eachindex(a);s+=a[i]*b[i];end;s)

function dirichlet_upper(N::Int; nt::Int=48, npanel::Int=48, nodes::Vector{IV}, wts::Vector{IV}, verbose=true)
    T0=time()
    ta=@elapsed (H,D,nx)=assemble_dirichlet(N; nt=nt, npanel=npanel, nodes=nodes, wts=wts)
    verbose && @printf("  [dirichlet assemble N=%d] D=%d modes/axis=%d (%.1fs)\n",N,D,nx,ta); flush(stdout)
    Hm = mid.(H)
    es = eigen(Symmetric(Hm))
    v = es.vectors[:,1]; v ./= norm(v); vI=interval.(v)
    Hv = imatvec(H,vI); num=idot(vI,Hv); den=idot(vI,vI)
    RQ = num/den
    ub = sup(RQ)                            # rigorous upper bound on lambda1^D >= lambda1(R^3)
    verbose && @printf("  [dirichlet] float lambda1^D_N=%.12f\n", es.values[1])
    verbose && @printf("  [dirichlet] RQ_interval=[%.12f,%.12f] => UPPER=%.12f (den w=%.2e)\n",
                       inf(RQ),sup(RQ),ub,diam(den)); flush(stdout)
    return Dict{String,Any}("N"=>N,"D"=>D,"modes_per_axis"=>nx,
        "float_lambda1D_N"=>es.values[1],
        "RQ_interval"=>[inf(RQ),sup(RQ)],"dirichlet_upper_rigorous"=>ub,
        "RQ_width"=>diam(RQ),"mass_den"=>[inf(den),sup(den)],
        "t_assemble"=>ta,"wall_total"=>time()-T0)
end
end # module
