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
    spec = zeros(Float64,nη,nw)
    
    for i=1:nη
        for j=1:nw
            u = split(dataset[1+(i-1)*nη+j])
            num = parse(Float64,u[3])
            spec[i,j] = num
        end
    end

    return spec
end

U = 20.0
nη = 100
nw = 200
ηs = range(0.0,6.0,length=nη)
ws = range(0.1,4*U,length=nw)
spec = readdata()

fig = heatmap(ηs,ws,spec)
savefig(fig,"opcon_4site4photo.png")


