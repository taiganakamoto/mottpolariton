using LinearAlgebra
using Plots

function readdata()
    fn = "spec.txt"
    dataset = readlines(fn)
    numdata = countlines(fn)
    U = 20.0
    ηs = []
    ws = []
    spec_lan = zeros(Float64,nw,nη)
    
    for i=1:nη 
        for j=1:nw
            u = split(dataset[(i-1)*nw+j+1])
            spec_lan[j,i] = parse(Float64,u[3])
            if i == 1
                push!(ws,parse(Float64,u[2]))
            end
            if j == 1
                push!(ηs,parse(Float64,u[1]))
            end
        end
    end

    return ηs, ws, spec_lan
end

nsite = 5
nelec = nsite
ndown = div(nsite,2)
N_c = 4
U = 20.0
μ = U/2
ω = U
nη = 100
nw = 200

ηs, ws, spec_lan = readdata()

fig = heatmap(ηs,ws,spec_lan,c=:thermal,clims=(0,50),title="N_e=$nelec,N_c=$N_c,ω_c=$ω,U=$U",colorbar_title="ω*σ(ω)",xlabel="η",ylabel="ω")

savefig(fig,"opcon_$(nelec)e$(N_c)p$(U)U$(ω)w.png")