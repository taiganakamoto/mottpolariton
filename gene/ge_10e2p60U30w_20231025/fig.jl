using LinearAlgebra
using Plots

function readdata()
    fn = "gene.txt"
    dataset = readlines(fn)
    nη = countlines(fn) - 1
    ηs = zeros(Float64,nη)
    spec = zeros(Float64,nη)
    
    for i=1:nη 
        u = split(dataset[i+1])
        spec[i] = parse(Float64,u[2])
        ηs[i] = parse(Float64,u[1])
    end

    return ηs, spec
end

nsite = 10
nelec = nsite
ndown = div(nsite,2)
N_c = 2
U = parse(Float64,ARGS[1])
μ = U/2
ω = U/2
δ = 0.1

function ge_nocav(ηs,U,nelec)
    nη = length(ηs)
    u = zeros(Float64,nη)
    for i=1:nη
        u[i] = -4*log(2)*nelec/U - U*nelec/4
    end
    return u
end

function ge_cav(ηs,U,nelec)
    nη = length(ηs)
    u = zeros(Float64,nη)
    for i=1:nη
        u[i] = -4*(1+ηs[i]^2/nelec/sqrt(1+ηs[i]^2))*log(2)*nelec/U - U*nelec/4
    end
    return u
end

ηs, spec = readdata()

fig = plot!(ηs,spec,title="N_e=$nelec,N_c=$N_c,ω_c=$ω,U=$U",xlabel="η",ylabel="Ground energy",label="numerical plot")
fig = plot!(ηs,ge_nocav(ηs,U,nelec),label="Heisenberg gene")
fig = plot!(ηs,ge_cav(ηs,U,nelec),label="modified Heisenberg gene")

savefig(fig,"ge_$(nelec)e$(N_c)p$(U)U$(ω)w.png")