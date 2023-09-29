#試行ベクトルをランダムにしてLanczos法を実行
function lanczos(hamil, steps)
    Random.seed!(5)
    ndim = size(hamil,1)
    v = zeros(ComplexF64, ndim, steps+1)
    alpha = zeros(Float64, steps)
    beta = zeros(Float64, steps)
    v[:,1] = rand(ComplexF64,ndim)

    # Lanczos steps
    # first step
    beta[1] = 1.0
    v[:,1] /= norm(v[:,1])
    alpha[1] = real(v[:,1]' * hamil * v[:,1])
    v[:,2] = hamil * v[:,1] - alpha[1] * v[:,1]
    # jth steps
    for j in 2:steps
        beta[j] = norm(v[:,j])
        v[:,j] /= beta[j]
        alpha[j] = real(v[:,j]' * hamil * v[:,j])
        v[:,j+1] = hamil * v[:,j] - alpha[j] * v[:,j] - beta[j] * v[:,j-1]
    end
    T = SymTridiagonal(alpha,beta[2:steps])
    es, vs = eigen(T)
    vg = v[:,1:steps] * vs[:,1]
    return es[1], vg
end

#試行ベクトルをv1として、Lanczos法を実行
function lanczos(hamil, v1, steps)
    ndim = length(v1)
    v = zeros(ComplexF64, ndim, steps+1)
    alpha = zeros(Float64, steps)
    beta = zeros(Float64, steps)
    v[:, 1] = v1

    # Lanczos steps
    # first step
    beta[1] = 1.0
    v[:,1] /= norm(v[:,1])
    alpha[1] = real(v[:,1]' * hamil * v[:,1])
    v[:,2] = hamil * v[:,1] - alpha[1] * v[:,1]
    # jth steps
    for j in 2:steps
        beta[j] = norm(v[:,j])
        v[:,j] /= beta[j]
        alpha[j] = real(v[:,j]' * hamil * v[:,j])
        v[:,j+1] = hamil * v[:,j] - alpha[j] * v[:,j] - beta[j] * v[:,j-1]
    end
    return alpha, beta
end

function frac(n, z, alpha, beta)
    if n < length(alpha)
        return beta[n]^2/(z-alpha[n]-frac(n+1,z,alpha,beta))
    else
        return 0
    end
end

function opcon(hamil,steps,ws,δ,nsite,nelec,N_c)
    eg, vg = lanczos(hamil,steps)
    J = currents(nsite,nelec,N_c)
    v1 = J*vg
    nm = norm(v1)
    v1 /= nm
    alpha, beta = lanczos(hamil,v1,steps)
    nw = length(ws)
    spec = zeros(Float64,nw)

    for i=1:nw
        z = ws[i] + δ*1im + eg
        X = frac(1,z,alpha,beta)
        spec[i] = -nm^2*imag(X)
    end

    return spec
end