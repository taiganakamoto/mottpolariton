function i2sites(i,nsite)
    n = 2*nsite
    sites = zeros(Bool,n)
    ii = i
    for i=1:n
        sites[i] = ii % 2
        ii = div(ii-sites[i],2)
    end
    
    return sites
end

function Base.display(sites::Array{Bool,1})
    n = length(sites)
    print("|")
    for i=n:-1:1
        if i == div(n,2)
            print(";")
        end
        if sites[i]
            print("1")
        else
            print("0")
        end

    end
    println(">")
end

i = 4
println(i)
nsite = 4
a = i2sites(i,nsite)
display(a)
