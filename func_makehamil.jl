#状態を01で表現。左半分がダウンスピン、右半分がアップスピン
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
function combi(n,r)
    if (r == 0 || r == n)
        return 1
    else
        return combi(n-1,r-1) + combi(n-1,r)
    end
end

#粒子数を固定した状態と全状態を行き来する行列を生成
function n_full(nsite,nelec)
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

#c^{\dagger}_{ith,ispin}c_{jth,jspin} * |ip><jp|の行列表示
function make_cdc(ith,ispin,jth,jspin,ip,jp,nsite,nelec,N_c)
    isite = (ispin-1)*nsite + ith
    jsite = (jspin-1)*nsite + jth
    ntofull, fullton = n_full(nsite,nelec)
    ndim = length(ntofull)
    Ndim = ndim*(N_c+1)
    mat_cdc = spzeros(ComplexF64,Ndim,Ndim)
    
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
                mat_cdc[ip*ndim+m,jp*ndim+k] = sign1*sign2
            end
        end        
    end
    return mat_cdc
end

#N_cをフォトンのカットオフとして、全ハミルトニアン。フォトンのセクターが0から始まる事に注意。U=0は注意。
function make_hamil(μ,U,η,ω,nsite,nelec,N_c)
    ndim = combi(2*nsite,nelec)
    Ndim = ndim*(N_c+1)
    hamil = spzeros(ComplexF64,Ndim,Ndim)

    Ω = ω#*sqrt(1+η^2)
    ζ = η/sqrt(nelec)#/(1+η^2)^(1/4)
    
    for ip=0:N_c
        
        #フォトン数ipのセクターでの単位行列
        Uni = spzeros(ComplexF64,Ndim,Ndim)
        for i=1:ndim
            Uni[ip*ndim+i,ip*ndim+i] += complex(1.0,0.0)
        end

        #対角要素にフォトン数に応じたエネルギーを埋める。Lanczos法のために、0フォトンstateのエネルギーを絶対値最小にしておく
        hamil += (ip-N_c)*Ω*Uni

        #フォトン数が保存するセクターに行列要素を埋める
        jp = ip
        for ith=1:nsite
            jth = ith+1
            jth += ifelse(jth > nsite, -nsite, 0)
            for ispin=1:2
                jspin = ispin
                hamil += -make_cdc(ith,ispin,jth,jspin,ip,jp,nsite,nelec,N_c)
            end
            
            jth = ith-1
            jth += ifelse(jth < 1, nsite, 0)
            for ispin=1:2
                jspin = ispin
                hamil += -make_cdc(ith,ispin,jth,jspin,ip,jp,nsite,nelec,N_c)
            end

            hamil += U*(make_cdc(ith,2,ith,2,ip,jp,nsite,nelec,N_c)-(μ/U)*Uni)*(make_cdc(ith,1,ith,1,ip,jp,nsite,nelec,N_c)-(μ/U)*Uni)
        end

        #フォトン数が異なるセクターに行列要素を埋める
        jp = ip+1
        if jp <= N_c
            for ith=1:nsite
                jth = ith+1
                jth += ifelse(jth > nsite, -nsite, 0)
                for ispin=1:2
                    jspin = ispin
                    hamil += 1im*ζ*sqrt(jp)*make_cdc(ith,ispin,jth,jspin,ip,jp,nsite,nelec,N_c)
                end
                
                jth = ith-1
                jth += ifelse(jth < 1, nsite, 0)
                for ispin=1:2
                    jspin = ispin
                    hamil += -1im*ζ*sqrt(jp)*make_cdc(ith,ispin,jth,jspin,ip,jp,nsite,nelec,N_c)
                end
            end
        end

        jp = ip-1
        if jp >= 0
            for ith=1:nsite
                jth = ith+1
                jth += ifelse(jth > nsite, -nsite, 0)
                for ispin=1:2
                    jspin = ispin
                    hamil += 1im*ζ*sqrt(ip)*make_cdc(ith,ispin,jth,jspin,ip,jp,nsite,nelec,N_c)
                end
                
                jth = ith-1
                jth += ifelse(jth < 1, nsite, 0)
                for ispin=1:2
                    jspin = ispin
                    hamil += -1im*ζ*sqrt(ip)*make_cdc(ith,ispin,jth,jspin,ip,jp,nsite,nelec,N_c)
                end
            end
        end
    end
    return hamil
end

#dimension less current
function currents(nsite,nelec,N_c)
    ndim = combi(2*nsite,nelec)
    Ndim = ndim*(N_c+1)
    currents = spzeros(ComplexF64,Ndim,Ndim)

    for ip=0:N_c
        jp = ip
        for ith=1:nsite
            jth = ith+1
            jth += ifelse(jth > nsite, -nsite, 0)
            for ispin=1:2
                jspin = ispin
                currents += 1im*make_cdc(ith,ispin,jth,jspin,ip,jp,nsite,nelec,N_c)
            end
            
            jth = ith-1
            jth += ifelse(jth < 1, nsite, 0)
            for ispin=1:2
                jspin = ispin
                currents += -1im*make_cdc(ith,ispin,jth,jspin,ip,jp,nsite,nelec,N_c)
            end
        end
    end
    return currents
end