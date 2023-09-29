using LinearAlgebra
using Plots
using SparseArrays
using IterativeSolvers

include("func_hubbard.jl")

nsite = 4
nelec = nsite

nU = 100
Us = range(-10,length=nU,stop=10)
μs = Us/2
es = zeros(Float64,nU)
es2 = zeros(Float64,nU,70)
for i=1:nU
    hamil = make_hamil(μs[i],Us[i],nsite,nelec)
    hamiltonian = Array(hamil)
    e,v = eigen(hamiltonian)
    es[i] = e[2]-e[1]
    es2[i,:] = e[:]
    println(Us[i],"\t",es[i])
end
fig = plot(Us,es2,labels="",xlabel="U",ylabel="all energy")
savefig(fig,"testfig/Fig5_me.png")