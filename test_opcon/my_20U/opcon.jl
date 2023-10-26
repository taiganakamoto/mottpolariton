using LinearAlgebra
using SparseArrays
using Random

include("func_makehamil2.jl")
include("func_lanczos2.jl")

nsite = 8
nelec = nsite
N_c = 0
U = parse(Float64,ARGS[1])
μ = U/2
ω = U/2
nw = 200
δ = 0.1

ndim = combi(2*nsite,nelec)
Ndim = ndim*(N_c+1)
η = 0
ws = Vector(0:0.2:39.8)
spec = zeros(Float64,nw)

hamil = make_hamil(μ,U,η,ω,nsite,nelec,N_c)
spec = opcon(hamil,150,ws,δ,nsite,nelec,N_c)

out = open("spec.txt","w")
println(out,"η","\t","ω","\t","spec_lanczos")

for k=1:nw
    println(out,ws[k],"\t",spec[k])
end

close(out)