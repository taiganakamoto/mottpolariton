using LinearAlgebra
using Plots
using SparseArrays
using Random

include("func_makehamil.jl")
include("Lanczos.jl")

function combi(n,r)
    if (r == 0 || r == n)
        return 1
    else
        return combi(n-1,r-1) + combi(n-1,r)
    end
end

#
nsite = 4
nelec = nsite
N_c = 4
ndim = combi(2*nsite,nelec)
Ndim = ndim*(N_c+1)

U = 20.0
μ = U/2
nη = 100
ηs = range(0.0,6.0,length=nη) 
ω = U
ess = zeros(Float64,nη,Ndim)

for i=1:nη
    hamil = make_hamil(μ,U,ηs[i],ω,nsite,nelec,N_c)
    eg, vg = lanczos(hamil,50)
    hamil = Array(hamil)
    es,vs = eigen(hamil)
    ess[i,:] = real(es[:]) 
    println(real(es[1]),"\t",eg)
end
fig = plot(ηs,ess,title="N_c=4,ω=U=20",labels="",xlabel="η",ylabel="all energy")
savefig(fig,"Fig_all.png")
#

#=
nsite = 4
nelec = nsite
N_c = 4
ndim = combi(2*nsite,nelec)
Ndim = ndim*(N_c+1)

U = 20.0
μ = U/2
nη = 100
ηs = range(0.0,6.0,length=nη) 
ω = U
for i=1:nη
    hamil = make_hamil(μ,U,ηs[i],ω,nsite,nelec,N_c)
    eg, vg, alpha, beta = lanczos(hamil,200)
    println(real(eg))
end
=#