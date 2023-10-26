using LinearAlgebra
using SparseArrays
using Random
using Dates
using Base.Threads

include("func_makehamil.jl")
include("func_lanczos.jl")

nsite = parse(Int64,ARGS[1])
nelec = nsite
ndown = div(nsite,2)
N_c = parse(Int64,ARGS[2])
U = 20.0
μ = U/2
ω = U/2
nη = 28*4
nw = 200
δ = 0.05

ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
Ndim = ndim*(N_c+1)
ηs = range(0.0,4.0,length=nη)
ws = range(0.1,2*U,length=nw)
spec = zeros(Float64,nw,nη)

Threads.@threads for i=1:nη
    hamil = make_hamil(μ,U,ηs[i],ω,nsite,nelec,ndown,N_c)
    spec[:,i] = opcon(hamil,150,ws,δ,nsite,nelec,ndown,N_c)
end

out = open("spec.txt","w")
println(out,"η","\t","ω","\t","spec_lanczos")

for i=1:nη
    for k=1:nw
        println(out,ηs[i],"\t",ws[k],"\t",spec[k,i])
    end
end

close(out)

