using LinearAlgebra
using SparseArrays
using Random
using Base.Threads

include("func_makehamil2.jl")
include("func_lanczos2.jl")

nsite = parse(Int64,ARGS[1])
nelec = nsite
ndown = div(nsite,2)
N_c = 2
U = 20.0
μ = U/2
ω = U/2
nη = 28*4
nw = 200
δ = 0.05

ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
Ndim = ndim*(N_c+1)
ηs = range(0.0,U,length=nη)
ges = zeros(Float64,nη)

Threads.@threads for i=1:nη
    hamil = make_hamil(μ,U,ηs[i],ω,nsite,nelec,ndown,N_c)
    eg, v = lanczos(hamil,200)
    Ω = ω*sqrt(1+ηs[i]^2)
    ges[i] = eg + N_c*Ω
end

out = open("gene.txt","w")
println(out,"η","\t","gene")

for i=1:nη
    println(out,ηs[i],"\t",ges[i])
end