using LinearAlgebra
using Plots

include("func_nagai.jl")

nsite = 4
nelec = nsite

nU = 100
Us = range(-10,length=nU,stop=10)
μs = Us/2
es = zeros(Float64,nU)
es2 = zeros(Float64,nU,70)
for i=1:nU
    hamiltonian = make_hamiltonian(μs[i],Us[i],nsite)
    hamiltonian =reduce_hamiltonian(hamiltonian,nsite,nelec)
    e,v = eigen(hamiltonian)
    es[i] = e[2]-e[1]
    es2[i,:] = e[:] .+ μs[i]^2*nelec/Us[i]
    println(Us[i],"\t",es[i])
end
fig = plot(Us,es2,labels="",xlabel="U",ylabel="all energy")
savefig(fig,"testfig/Fig5_nagai.png")