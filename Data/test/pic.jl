using LinearAlgebra
using Plots

function readdata()
    fn = "spec.txt"
    dataset = readlines(fn)
    numdata = countlines(fn)
    U = 20.0
    nη = 100
    nw = 200
    ηs = range(0.0,6.0,length=nη)
    ws = range(0.1,4*U,length=nw)
    spec_exa = zeros(Float64,nw,nη)
    spec_lan = zeros(Float64,nw,nη)
    
    for i=1:nη 
        for j=1:nw
            u = split(dataset[(i-1)*nw+j+1])
            spec_exa[j,i] = parse(Float64,u[3])
            spec_lan[j,i] = parse(Float64,u[4])
        end
    end

    return spec_exa, spec_lan
end

nsite = 4
nelec = nsite
N_c = 3
U = 20.0
nη = 100
nw = 200
ηs = range(0.0,6.0,length=nη)
ws = range(0.1,4*U,length=nw)
spec_exa,spec_lan = readdata()

fig = heatmap(ηs,ws,spec_exa,c=:thermal,title="N_e=$nelec,N_c=$N_c,ω_c=U=$U",colorbar_title="ω*σ(ω)",xlabel="η",ylabel="ω")
fig2 = heatmap(ηs,ws,spec_lan,c=:thermal,title="N_e=$nelec,N_c=$N_c,ω_c=U=$U",colorbar_title="ω*σ(ω)",xlabel="η",ylabel="ω")

savefig(fig,"opcon_4e3p_exact.png")
savefig(fig2,"opcon_4e3p_lanczos.png")