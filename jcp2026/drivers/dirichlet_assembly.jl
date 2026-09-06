# dirichlet_assembly.jl — sine-basis Dirichlet interval assembly for H2+ on Omega2,
# generalized to arbitrary sine-parity symmetry sectors (for locating & certifying
# the 2nd Dirichlet eigenvalue as a Lehmann-Goerisch separator).
#
# Sine basis psi_n(x) = (1/sqrt(L)) sin(n*pi*(x+L)/(2L)), n=1,2,...  (vanish at x=+-L).
#   kinetic eigenvalue (n*pi/(2L))^2 ; L2-orthonormal (mass = identity).
#   Coulomb sine moment 0.5*(1/L)*(G[|n-m|] - G[n+m])  (product-to-sum: MINUS).
# Parity in x about 0: sin(n*pi*(x+L)/(2L)) is EVEN-in-x for ODD n, ODD-in-x for EVEN n.
#   x-sector :odd  -> odd sine modes n=1,3,5,... (even-in-x)  => Dirichlet GROUND sector
#   x-sector :even -> even sine modes n=2,4,6,... (odd-in-x)
# y,z always :odd (even-in-y/z) for the ground; a sector flips one axis to reach lambda2^D.
module DirichletAssembly
using LinearAlgebra, Printf
using IntervalArithmetic: Interval, interval, mid, inf, sup, diam
include("moments_verified.jl")
using .MomentsVerified
const IV = Interval{Float64}
const LX, LY, LZ = 20.0, 16.0, 16.0
const AX = 2.0
t_star_() = 1.0

build_t_grid(nt::Int) = begin
    k=1:(nt-1); β=collect(k ./ sqrt.(4.0 .*k.^2 .-1.0))
    J=SymTridiagonal(zeros(nt),β); vals,vecs=eigen(J)
    x=vals; w=2.0 .*(vecs[1,:].^2); u=0.5 .*(x .+1); wu=0.5 .*w
    (u ./(1 .-u), wu ./(1 .-u).^2)
end

function sine_moment_matrix_1d(mmax::Int, L::Float64, s::Float64, t::Float64,
                               nodes::Vector{IV}, wts::Vector{IV}; npanel::Int=96, t_star::Float64=1.0)
    fmax = 2*mmax
    G = Vector{IV}(undef, fmax+1)
    for f in 0:fmax
        κ = f*pi/(2L)
        G[f+1] = MomentsVerified.moment_verified(κ, t, L, s; t_star=t_star, nodes=nodes, wts=wts, npanel=npanel)
    end
    Ninv = interval(1.0)/interval(L); half = interval(0.5)
    B = Matrix{IV}(undef, mmax, mmax)
    @inbounds for n in 1:mmax, m in 1:mmax
        B[n,m] = half*Ninv*(G[abs(n-m)+1] - G[n+m+1])
    end
    return B
end

# sine modes of a given parity among 1..mmax
sector_modes(mmax::Int, par::Symbol) = par===:odd ? collect(1:2:mmax) : collect(2:2:mmax)

# assemble Dirichlet POTENTIAL matrix P and KINETIC diagonal (separately) for sector (px,py,pz).
function assemble_PK(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                  nt::Int=48, npanel::Int=96, nodes::Vector{IV}, wts::Vector{IV})
    mmax = N+1
    t, wt = build_t_grid(nt)
    Ix = sector_modes(mmax,px); Iy = sector_modes(mmax,py); Iz = sector_modes(mmax,pz)
    nx=length(Ix); ny=length(Iy); nz=length(Iz); D=nx*ny*nz
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
    νx=(interval(pi)/(interval(2.0)*interval(LX)))^2
    νy=(interval(pi)/(interval(2.0)*interval(LY)))^2
    νz=(interval(pi)/(interval(2.0)*interval(LZ)))^2
    kin = Vector{IV}(undef, D)
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        idx=(cz-1)*ny*nx+(cy-1)*nx+cx
        kin[idx] = νx*interval(Ix[cx]^2)+νy*interval(Iy[cy]^2)+νz*interval(Iz[cz]^2)
    end
    return P, kin, D, (nx,ny,nz)
end

# form full interval Hamiltonian H = P + diag(kin_scale*kin) + sigma*I
function form_H(P::Matrix{IV}, kin::Vector{IV}; kin_scale::IV=interval(1.0), sigma::IV=interval(0.0))
    D=length(kin); H=copy(P)
    @inbounds for i in 1:D
        H[i,i] = H[i,i] + kin_scale*kin[i] + sigma
    end
    return H
end

# convenience: assemble full H for a sector (kin_scale=1)
function assemble(N::Int, px::Symbol, py::Symbol, pz::Symbol;
                  nt::Int=48, npanel::Int=96, nodes::Vector{IV}, wts::Vector{IV})
    P,kin,D,dims = assemble_PK(N,px,py,pz; nt=nt,npanel=npanel,nodes=nodes,wts=wts)
    return form_H(P,kin), D, dims
end

imatvec(M::Matrix{IV}, x::Vector{IV}) = begin
    D=length(x); y=Vector{IV}(undef,D)
    @inbounds for i in 1:D
        acc=interval(0.0); @simd for j in 1:D; acc+=M[i,j]*x[j]; end; y[i]=acc
    end; y
end
idot(a::Vector{IV},b::Vector{IV})=(s=interval(0.0);@inbounds for i in eachindex(a);s+=a[i]*b[i];end;s)
end # module
