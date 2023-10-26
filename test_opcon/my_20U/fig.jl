using LinearAlgebra
using Plots

function readdata1()
    fn = "spec.txt"
    dataset = readlines(fn)
    nw = countlines(fn) - 1
    ws = zeros(Float64,nw)
    spec = zeros(Float64,nw)
    
    for i=1:nw
        u = split(dataset[i+1])
        ws[i] = parse(Float64,u[1])
        spec[i] = parse(Float64,u[2])
    end

    return ws, spec
end

function readdata2()
    fn = "zvo_DynamicalGreen.dat"
    dataset = readlines(fn)
    nw = countlines(fn)
    ws = zeros(Float64,nw)
    spec = zeros(Float64,nw)
    
    for i=1:nw
        u = split(dataset[i])
        ws[i] = parse(Float64,u[1])
        spec[i] = -parse(Float64,u[4])
    end

    return ws, spec
end

nelec = 8
U = parse(Float64,ARGS[1])
δ = 0.1

ws1, spec1 = readdata1()
ws2, spec2 = readdata2()

fig = plot!(ws1,spec1,title="N_e=$nelec,U=$U",xlabel="ω",ylabel="ω*σ(ω)",markershape=:circle,mc="blue",label="my code")
fig = plot!(ws2,spec2,title="N_e=$nelec,U=$U",xlabel="ω",ylabel="ω*σ(ω)",markershape=:circle,mc="red",label="H φ")

savefig(fig,"hphi&my_$(nelec)e$(U)U.png")