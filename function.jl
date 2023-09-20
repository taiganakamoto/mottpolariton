#状態を01で表現。左半分がダウンスピン、右半分がアップスピン
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

#状態をケット表示して表す関数
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

#どのサイトがどの状態の番号なのかを調べる関数
function sites2i(sites)
    n = length(sites)
    ii = 0
    for i=1:n
        ii += sites[i]*2^(i-1)
    end
    return ii
end

#サイトに電子がいるかどうかの関数
function check_site(isite,sites)
    return sites[isite]
end

#考えているサイトより左の粒子の数に応じて符号を与える関数
function check_sign(isite,sites)
    n = length(sites)
    sign = 1
    for jsite = n:-1:isite+1
        if sites[jsite]
            sign *= -1
        end
    end
    return sign
end

#i番目の消滅演算子
function make_c(ith,ispin,nsite)
    isite = (ispin-1)*nsite +ith
    ndim = 4^nsite
    mat_c = zeros(Int,ndim,ndim)
    
    for j=1:ndim
        sites = i2sites(j,nsite)
        if check_site(isite,sites)
            sign = check_sign(isite,sites)
            sites[isite] = false
            i = sites2i(sites)
            mat_c[i+1,j+1] = sign
        end        
    end
    return mat_c
end

#生成消滅演算子の組
function make_cs(nsite)
    vec_c = []
    vec_cdag = []
    for ispin = 1:2
        for ith=1:nsite        
            mat_c = make_c(ith,ispin,nsite)
            push!(vec_c,mat_c)
            push!(vec_cdag,mat_c')
        end
    end
    
    return vec_c,vec_cdag
end

#fullのハミルトニアンの構成
function make_hamiltonian(μ,U,nsite)
    ndim = 4^nsite
    hamiltonian = zeros(Float64,ndim,ndim)
    vec_c,vec_cdag = make_cs(nsite)
    for i=1:nsite
        j = i+1
        j += ifelse(j > nsite,-nsite,0)
        
        if 1 <= j <= nsite
        for ispin=1:2
            ii = (ispin-1)*nsite + i
            jj = (ispin-1)*nsite + j
            hamiltonian += - vec_cdag[ii]*vec_c[jj]
        end
        end
        
        j = i-1
        j += ifelse(j < 1,nsite,0)
        if 1 <= j <= nsite
        for ispin=1:2
            ii = (ispin-1)*nsite + i
            jj = (ispin-1)*nsite + j
            hamiltonian += - vec_cdag[ii]*vec_c[jj]
        end
        end

    
        j = i
        for ispin=1:2
            ii = (ispin-1)*nsite + i
            jj = (ispin-1)*nsite + j
            hamiltonian += -μ*vec_cdag[ii]*vec_c[jj]
        end
    
        iup = i
        idown = nsite + i
        hamiltonian += U*vec_cdag[iup]*vec_c[iup]*vec_cdag[idown]*vec_c[idown]
   
    
    
    end
    return hamiltonian
end

function make_mapping(nsite,nelec)
    imap = zeros(Int,4^nsite)
    icount = 0
    for i=1:4^nsite
        sites = i2sites(i-1,nsite)
        if sum(sites) == nelec
            icount += 1
            imap[i] = icount
        end
    end
    return imap,icount
end

function reduce_hamiltonian(ham,nsite,nelec)
    ndim = size(ham)[1]     
    @assert ndim == 4^nsite
    
    imap,numdim = make_mapping(nsite,nelec)
    ham_reduce = zeros(Float64,numdim,numdim)
    
    for i=1:ndim
        ip = imap[i]
        
        if ip  != 0
            for j=1:ndim
                jp = imap[j] 
                if jp != 0
                    ham_reduce[ip,jp] = ham[i,j]
                end            
            end
        end
    end
    return ham_reduce
    
end