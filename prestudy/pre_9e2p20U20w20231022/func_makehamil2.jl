#状態を01で表現。上半分がダウンスピン、下半分がアップスピン
@inline function i2sites(i,nsite)
    n = 2*nsite
    sites = zeros(Bool,n)
    ii = i
    for i=1:n
        sites[i] = ii % 2
        ii = div(ii-sites[i],2)
    end  
    return sites
end

#どのサイトがどの状態の番号なのかを調べる関数
@inline function sites2i(sites)
    n = length(sites)
    ii = 0
    for i=1:n
        ii += sites[i]*2^(i-1)
    end
    return ii
end

#サイトに電子がいるかどうかの関数
@inline function check_site(isite,sites)
    return sites[isite]
end

#考えているサイトより左の粒子の数に応じて符号を与える関数
@inline function check_sign(isite,sites)
    n = length(sites)
    sign = 1
    for jsite = n:-1:isite+1
        if sites[jsite]
            sign *= -1
        end
    end
    return sign
end

#組み合わせ計算
@inline function combi(n,r)
    if (r == 0 || r == n)
        return 1
    else
        return combi(n-1,r-1) + combi(n-1,r)
    end
end

#粒子数を固定した状態と全状態を行き来する行列を生成
@inline function n_full(nsite,nelec)
    count = 1
    ndim = combi(2*nsite,nelec)
    fulldim = 4^nsite
    ntofull = zeros(Int64,ndim)
    fullton = zeros(Int64,fulldim)
    for i = 1:4^nsite
        if sum(i2sites(i,nsite)) == nelec
            ntofull[count] = i
            fullton[i] = count
            count += 1
        end
    end
    return ntofull, fullton
end

@inline function n_full(nsite,nelec,ndown)
    count = 1
    ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
    fulldim = 4^nsite
    ntofull = zeros(Int64,ndim)
    fullton = zeros(Int64,fulldim)
    for i = 1:fulldim
        v = i2sites(i,nsite)
        len_dw = div(length(v),2)
        if sum(v) == nelec && sum(v[1:len_dw]) == ndown
            ntofull[count] = i
            fullton[i] = count
            count += 1
        end
    end
    return ntofull, fullton
end

#c^{\dagger}_{ith,ispin}c_{jth,jspin}の行列表示
@inline function make_cdc(ith,ispin,jth,jspin,nsite,nelec)
    isite = (ispin-1)*nsite + ith
    jsite = (jspin-1)*nsite + jth
    ntofull, fullton = n_full(nsite,nelec)
    ndim = length(ntofull)
    mat_cdc = spzeros(Float64,ndim,ndim)
    
    for k=1:ndim
        k_full = ntofull[k]
        sites = i2sites(k_full,nsite)
        if check_site(jsite,sites)
            sign1 = check_sign(jsite,sites)
            sites[jsite] = false
            if ! check_site(isite,sites)
                sign2 = check_sign(isite,sites)
                sites[isite] = true
                m_full = sites2i(sites)
                m = fullton[m_full]
                mat_cdc[m,k] = sign1*sign2
            end
        end        
    end
    return mat_cdc
end

@inline function make_cdc(ith,ispin,jth,jspin,nsite,nelec,ndown)
    isite = (ispin-1)*nsite + ith
    jsite = (jspin-1)*nsite + jth
    ntofull, fullton = n_full(nsite,nelec,ndown)
    ndim = length(ntofull)
    mat_cdc = spzeros(Float64,ndim,ndim)
    
    for k=1:ndim
        k_full = ntofull[k]
        sites = i2sites(k_full,nsite)
        if check_site(jsite,sites)
            sign1 = check_sign(jsite,sites)
            sites[jsite] = false
            if ! check_site(isite,sites)
                sign2 = check_sign(isite,sites)
                sites[isite] = true
                m_full = sites2i(sites)
                m = fullton[m_full]
                mat_cdc[m,k] = sign1*sign2
            end
        end        
    end
    return mat_cdc
end

#diagonal matrices of ones
function make_Uni(nsite,nelec)
    ndim = combi(2*nsite,nelec)
    Uni = spzeros(Float64,ndim,ndim)

    for i=1:ndim
        Uni[i,i] = 1.0
    end
    return Uni 
end

function make_Uni(nsite,nelec,ndown)
    ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
    Uni = spzeros(Float64,ndim,ndim)

    for i=1:ndim
        Uni[i,i] = 1.0
    end
    return Uni 
end

#dimension less current. In order no to use complex number, difinition is slightly changed
function make_current(nsite,nelec)
    ndim = combi(2*nsite,nelec)
    J = spzeros(Float64,ndim,ndim)

    for ith=1:nsite
        jth = ith+1
        jth += ifelse(jth > nsite, -nsite, 0)
        for ispin=1:2
            jspin = ispin
            J += make_cdc(ith,ispin,jth,jspin,nsite,nelec)
        end
            
        jth = ith-1
        jth += ifelse(jth < 1, nsite, 0)
        for ispin=1:2
            jspin = ispin
            J += -make_cdc(ith,ispin,jth,jspin,nsite,nelec)
        end
    end
    return J
end

function make_current(nsite,nelec,ndown)
    ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
    J = spzeros(Float64,ndim,ndim)

    for ith=1:nsite
        jth = ith+1
        jth += ifelse(jth > nsite, -nsite, 0)
        for ispin=1:2
            jspin = ispin
            J += make_cdc(ith,ispin,jth,jspin,nsite,nelec,ndown)
        end
            
        jth = ith-1
        jth += ifelse(jth < 1, nsite, 0)
        for ispin=1:2
            jspin = ispin
            J += -make_cdc(ith,ispin,jth,jspin,nsite,nelec,ndown)
        end
    end
    return J
end

#dimension less current including photon
function make_Current(nsite,nelec,N_c)
    ndim = combi(2*nsite,nelec)
    Ndim = ndim*(N_c+1)
    Current = spzeros(Float64,Ndim,Ndim)
    J = make_current(nsite,nelec)

    for ip=0:N_c
        jp = ip
        Current[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += J
    end
    return Current
end

function make_Current(nsite,nelec,ndown,N_c)
    ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
    Ndim = ndim*(N_c+1)
    Current = spzeros(Float64,Ndim,Ndim)
    J = make_current(nsite,nelec,ndown)

    for ip=0:N_c
        jp = ip
        Current[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += J
    end
    return Current
end

#hopping term
function make_hop(nsite,nelec)
    ndim = combi(2*nsite,nelec)
    T = spzeros(Float64,ndim,ndim)

    for ith=1:nsite
        jth = ith+1
        jth += ifelse(jth > nsite, -nsite, 0)
        for ispin=1:2
            jspin = ispin
            T += -make_cdc(ith,ispin,jth,jspin,nsite,nelec)
        end
            
        jth = ith-1
        jth += ifelse(jth < 1, nsite, 0)
        for ispin=1:2
            jspin = ispin
            T += -make_cdc(ith,ispin,jth,jspin,nsite,nelec)
        end
    end
    return T
end

function make_hop(nsite,nelec,ndown)
    ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
    T = spzeros(Float64,ndim,ndim)

    for ith=1:nsite
        jth = ith+1
        jth += ifelse(jth > nsite, -nsite, 0)
        for ispin=1:2
            jspin = ispin
            T += -make_cdc(ith,ispin,jth,jspin,nsite,nelec,ndown)
        end
            
        jth = ith-1
        jth += ifelse(jth < 1, nsite, 0)
        for ispin=1:2
            jspin = ispin
            T += -make_cdc(ith,ispin,jth,jspin,nsite,nelec,ndown)
        end
    end
    return T
end

#∑_{i}(n_{i,up}-μ/U)*(n_{i,down}-μ/U)
function make_nn(nsite,nelec,U,μ)
    ndim = combi(2*nsite,nelec)
    nn = spzeros(Float64,ndim,ndim)
    Uni = make_Uni(nsite,nelec)

    for ith=1:nsite
        nn += (make_cdc(ith,2,ith,2,nsite,nelec)-(μ/U)*Uni)*(make_cdc(ith,1,ith,1,nsite,nelec)-(μ/U)*Uni)
    end
    return nn
end

function make_nn(nsite,nelec,ndown,U,μ)
    ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
    nn = spzeros(Float64,ndim,ndim)
    Uni = make_Uni(nsite,nelec,ndown)

    for ith=1:nsite
        nn += (make_cdc(ith,2,ith,2,nsite,nelec,ndown)-(μ/U)*Uni)*(make_cdc(ith,1,ith,1,nsite,nelec,ndown)-(μ/U)*Uni)
    end
    return nn
end

#N_cをフォトンのカットオフとして、全ハミルトニアン。フォトンのセクターが0から始まる事に注意。U=0は注意。
function make_hamil(μ,U,η,ω,nsite,nelec,N_c)
    ndim = combi(2*nsite,nelec)
    Ndim = ndim*(N_c+1)
    hamil = spzeros(Float64,Ndim,Ndim)

    Ω = ω#*sqrt(1+η^2)
    ζ = η/sqrt(nelec)#/(1+η^2)^(1/4)

    Uni = make_Uni(nsite,nelec)
    J = make_current(nsite,nelec)
    T = make_hop(nsite,nelec)
    nn = make_nn(nsite,nelec,U,μ)
    
    for ip=0:N_c
        
        jp = ip
        #対角要素にフォトン数に応じたエネルギーを埋める。Lanczos法のために、0フォトンstateのエネルギーを絶対値最小にしておく。基底eneを見る際は、+N_c*Ωを全体に付ける必要がある。
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += (ip-N_c)*Ω*Uni
        
        #フォトン数が保存するセクターに行列要素を埋める。ホッピングと相互作用項
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += T
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += U*nn

        #フォトン数が異なるセクターに行列要素を埋める
        jp = ip+1
        if jp <= N_c
            hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] = -(ζ*sqrt(jp))*J
        end
        jp = ip-1
        if jp >= 0
            hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] = (ζ*sqrt(ip))*J
        end

    end
    return hamil
end

#N_cをフォトンのカットオフ、ndownをダウンスピンの数として、全ハミルトニアン。フォトンのセクターが0から始まる事に注意。U=0は注意。
function make_hamil(μ,U,η,ω,nsite,nelec,ndown,N_c)
    ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
    Ndim = ndim*(N_c+1)
    hamil = spzeros(Float64,Ndim,Ndim)

    Ω = ω#*sqrt(1+η^2)
    ζ = η/sqrt(nelec)#/(1+η^2)^(1/4)

    Uni = make_Uni(nsite,nelec,ndown)
    J = make_current(nsite,nelec,ndown)
    T = make_hop(nsite,nelec,ndown)
    nn = make_nn(nsite,nelec,ndown,U,μ)
    
    for ip=0:N_c
        
        jp = ip
        #対角要素にフォトン数に応じたエネルギーを埋める。Lanczos法のために、0フォトンstateのエネルギーを絶対値最小にしておく。基底eneを見る際は、+N_c*Ωを全体に付ける必要がある。
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += (ip-N_c)*Ω*Uni
        
        #フォトン数が保存するセクターに行列要素を埋める。ホッピングと相互作用項
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += T
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += U*nn

        #フォトン数が異なるセクターに行列要素を埋める
        jp = ip+1
        if jp <= N_c
            hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] = -(ζ*sqrt(jp))*J
        end
        jp = ip-1
        if jp >= 0
            hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] = (ζ*sqrt(ip))*J
        end

    end
    return hamil
end

function make_effhamil(μ,U,η,ω,nsite,nelec,N_c)
    ndim = combi(2*nsite,nelec)
    Ndim = ndim*(N_c+1)
    hamil = spzeros(Float64,Ndim,Ndim)

    Ω = ω*sqrt(1+η^2)
    ζ = η/sqrt(nelec)/(1+η^2)^(1/4)

    Uni = make_Uni(nsite,nelec)
    J = make_current(nsite,nelec)
    T = make_hop(nsite,nelec)
    nn = make_nn(nsite,nelec,U,μ)
    
    for ip=0:N_c
        
        jp = ip
        #対角要素にフォトン数に応じたエネルギーを埋める。Lanczos法のために、0フォトンstateのエネルギーを絶対値最小にしておく。基底状態からの差しか見ないなら、この取り扱いで問題ない。
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += (ip-N_c)*Ω*Uni
        
        #フォトン数が保存するセクターに行列要素を埋める。ホッピングと相互作用項
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += T
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += U*nn
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += (ζ^2/Ω)*J*J

        #フォトン数が異なるセクターに行列要素を埋める
        jp = ip+1
        if jp <= N_c
            hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] = -(ζ*sqrt(jp)/Ω)*(J*nn-nn*J)
        end
        jp = ip-1
        if jp >= 0
            hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] = (ζ*sqrt(ip)/Ω)*(J*nn-nn*J)
        end

    end
    return hamil
end

function make_effhamil(μ,U,η,ω,nsite,nelec,ndown,N_c)
    ndim = combi(nsite,ndown)*combi(nsite,nelec-ndown)
    Ndim = ndim*(N_c+1)
    hamil = spzeros(Float64,Ndim,Ndim)

    Ω = ω*sqrt(1+η^2)
    ζ = η/sqrt(nelec)/(1+η^2)^(1/4)

    Uni = make_Uni(nsite,nelec,ndown)
    J = make_current(nsite,nelec,ndown)
    T = make_hop(nsite,nelec,ndown)
    nn = make_nn(nsite,nelec,ndown,U,μ)
    
    for ip=0:N_c
        
        jp = ip
        #対角要素にフォトン数に応じたエネルギーを埋める。Lanczos法のために、0フォトンstateのエネルギーを絶対値最小にしておく。基底状態からの差しか見ないなら、この取り扱いで問題ない。
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += (ip-N_c)*Ω*Uni
        
        #フォトン数が保存するセクターに行列要素を埋める。ホッピングと相互作用項
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += T
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += U*nn
        hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] += (ζ^2/Ω)*J*J

        #フォトン数が異なるセクターに行列要素を埋める
        jp = ip+1
        if jp <= N_c
            hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] = -(ζ*sqrt(jp)/Ω)*(J*nn-nn*J)
        end
        jp = ip-1
        if jp >= 0
            hamil[(ip*ndim+1):(ip*ndim+ndim),(jp*ndim+1):(jp*ndim+ndim)] = (ζ*sqrt(ip)/Ω)*(J*nn-nn*J)
        end

    end
    return hamil
end