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

#c^{\dagger}_{ith,ispin}c_{jth,jspin}の行列表示
function make_cdc(ith,ispin,jth,jspin,nsite,nelec)
    isite = (ispin-1)*nsite + ith
    jsite = (jspin-1)*nsite + jth
    ntofull, fullton = n_full(nsite,nelec)
    ndim = length(ntofull)
    mat_cdc = spzeros(ndim,ndim)

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
                mat_cdc[m,k] = sign1 * sign2
            end
        end        
    end
    return mat_cdc
end

#ハミルトニアンの作成。ホッピングt=1。周期境界条件。
function make_hamil(μ,U,nsite,nelec)
    ndim = combi(2*nsite,nelec)
    hamil = spzeros(ndim,ndim)
    for ith=1:nsite
        jth = ith+1
        jth += ifelse(jth > nsite, -nsite, 0)
        for ispin=1:2
            jspin = ispin
            hamil += -make_cdc(ith,ispin,jth,jspin,nsite,nelec)
        end
        
        jth = ith-1
        jth += ifelse(jth < 1, nsite, 0)
        for ispin=1:2
            jspin = ispin
            hamil += -make_cdc(ith,ispin,jth,jspin,nsite,nelec)
        end

        hamil += U*(make_cdc(ith,2,ith,2,nsite,nelec)-(μ/U)*sparse(I,ndim,ndim))*(make_cdc(ith,1,ith,1,nsite,nelec)-(μ/U)*sparse(I,ndim,ndim))
    end
    
    return hamil
end


