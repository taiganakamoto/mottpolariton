using LinearAlgebra
using Plots

function readdata()
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

nelec = 8
U = parse(Float64,ARGS[1])
δ = 0.1

ws, spec = readdata()

fig = plot(ws,spec,title="N_e=$nelec,U=$U",xlabel="ω",ylabel="ω*σ(ω)")

savefig(fig,"my_$(nelec)e$(U)U.png")