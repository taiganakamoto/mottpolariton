using LinearAlgebra
using SparseArrays
using Random
using Dates

include("func_makehamil.jl")
include("func_lanczos.jl")

function combi(n,r)
    if (r == 0 || r == n)
        return 1
    else
        return combi(n-1,r-1) + combi(n-1,r)
    end
end

nsite = 6
nelec = nsite
N_c = 4
ndim = combi(2*nsite,nelec)
Ndim = ndim*(N_c+1)

U = 20.0
μ = U/2
ω = U
nη = 100
nw = 200
ηs = range(0.0,6.0,length=nη)
ws = range(0.1,2*U,length=nw)
δ = 0.05
spec = zeros(Float64,nw,nη)

out = open("spec.txt","w")
println(out,"η","\t","ω","\t","spec_lanczos")

for i=1:nη
    hamil = make_hamil(μ,U,ηs[i],ω,nsite,nelec,N_c)
    spec[:,i] = opcon(hamil,150,ws,δ,nsite,nelec,N_c)
    for k=1:nw
        println(out,ηs[i],"\t",ws[k],"\t",spec[k,i])
    end
end

close(out)