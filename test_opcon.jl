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
N_c = 3
ndim = combi(2*nsite,nelec)

U = 20.0
μ = U/2
ω = U
nη = 100
nw = 200
ηs = range(0.0,6.0,length=nη)
ws = range(0.1,2*U,length=nw)
δ = 0.01
spec = zeros(Float64,nw,nη)
spec2 = zeros(Float64,nw,nη)
J = currents(nsite,nelec,N_c)

for i=1:nη
    hamil = make_hamil(μ,U,ηs[i],ω,nsite,nelec,N_c)
    spec[:,i] = opcon(hamil,150,ws,δ,nsite,nelec,N_c)
    hamil2 = Array(hamil)
    es, vs = eigen(hamil2)
    Jvg = J*vs[:,1]
    for k=1:nw
        for j=1:ndim
            spec2[k,i] += real((vs[:,j]' * Jvg)^2)*imag(1/(ws[k] + es[1] - es[j] + δ*1im)) 
        end
    end
end
fig = heatmap(ηs,ws,spec,c=:thermal,title="ω=U=$U",colorbar_title="ω*σ(ω)",xlabel="η",ylabel="ω")
fig2 = heatmap(ηs,ws,spec2,c=:thermal,title="ω=η=$U",colorbar_title="ω*σ(ω)",xlabel="η",ylabel="ω")
savefig(fig,"Fig1_opcon_lanczos.png")
savefig(fig2,"Fig1_opcon_exact.png")
#


#=
nsite = 5
nelec = nsite
N_c = 0
ndim = combi(2*nsite,nelec)

nU = 100
Us = range(0.1,20.0,length=nU)
μs = Us/2
ω = 0.0
η = 0.0
nw = 200
ws = range(-1,40.0,length=nw)
δ = 0.05
spec = zeros(Float64,nw,nU)
spec2 = zeros(Float64,nw,nU)
J = currents(nsite,nelec,N_c)

for i=1:nU
    hamil = make_hamil(μs[i],Us[i],η,ω,nsite,nelec,N_c)
    spec[:,i] = opcon(hamil,100,ws,δ,nsite,nelec,N_c)
    hamil2 = Array(hamil)
    es, vs = eigen(hamil2)
    Jvg = J*vs[:,1]
    for k=1:nw
        for j=1:ndim
            spec2[k,i] += real((vs[:,j]' * Jvg)^2)*imag(1/(ws[k] + es[1] - es[j] + δ*1im)) 
        end
    end
end
fig = heatmap(Us,ws,spec,c=:thermal,title="ω=η=0",colorbar_title="ω*σ(ω)",xlabel="U",ylabel="ω")
fig2 = heatmap(Us,ws,spec2,c=:thermal,title="ω=η=0",colorbar_title="ω*σ(ω)",xlabel="U",ylabel="ω")
savefig(fig,"Fig0_opcon_lanczos.png")
savefig(fig2,"Fig0_opcon_exact.png")
=#

#=
nsite = 8
nelec = nsite
N_c = 0
ndim = combi(2*nsite,nelec)

U = 4.0
μ = U/2
ω = 0.0
η = 0.0
nw = 200
ws = range(0.1,2*U,length=nw)
δ = 0.05

hamil = make_hamil(μ,U,η,ω,nsite,nelec,N_c)
spec = opcon(hamil,100,ws,δ,nsite,nelec,N_c)

fig = plot(ws,spec)
savefig(fig,"Figsub_opcon_lanczos.png")
=#