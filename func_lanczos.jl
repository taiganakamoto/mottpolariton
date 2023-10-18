#試行ベクトルをランダムにしてLanczos法を実行
function lanczos(hamil, steps)
    Random.seed!(5)
    Ndim = size(hamil,1)
    alpha = zeros(Float64, steps)
    beta = zeros(Float64, steps)
    dum = zeros(ComplexF64,Ndim)
    v1 = rand(ComplexF64,Ndim)
    v2 = zeros(ComplexF64,Ndim)
    w1 = zeros(ComplexF64,Ndim)
    w2 = zeros(ComplexF64,Ndim)

    #initial state
    v1 = v1/norm(v1)
    w1 = v1 #copy initial state for calculating grandstate vector
    beta[1] = 1.0

    # lanczos steps
    for j in 1:steps-1
        v2 = v2 + hamil*v1
        alpha[j] = real(v1'*v2)
        dum = v1
        v1 = v2 - alpha[j]*v1
        v2 = dum
        beta[j+1] = norm(v1)
        v2 = -beta[j+1]*v2
        v1 = v1/beta[j+1]
    end

    #final step
    v2 = v2 + hamil*v1
    alpha[steps] = real(v1'*v2)

    #calculate tridiagonal matrices
    T = SymTridiagonal(alpha,beta[2:steps])
    es, vs = eigen(T)

    #calculate grandstate vector
    vg = zeros(ComplexF64,Ndim)
    vg += vs[1,1]*w1
    for j in 1:steps-1
        w2 = w2 + hamil*w1
        dum = w1
        w1 = w2 - alpha[j]*w1
        w2 = -beta[j+1]*dum
        w1 = w1/norm(beta[j+1])
        vg += vs[j+1,1]*w1
    end

    return es[1], vg
end

#規格化された試行ベクトルをv1として、Lanczos法を実行
function lanczos(hamil, v0, steps)
    Ndim = size(hamil,1)
    alpha = zeros(Float64, steps)
    beta = zeros(Float64, steps)
    dum = zeros(ComplexF64,Ndim)
    v1 = zeros(ComplexF64,Ndim)
    v2 = zeros(ComplexF64,Ndim)

    #initial state
    v1 = v0
    beta[1] = 1.0

    # lanczos steps
    for j in 1:steps-1
        v2 = v2 + hamil*v1
        alpha[j] = real(v1'*v2)
        dum = v1
        v1 = v2 - alpha[j]*v1
        v2 = dum
        beta[j+1] = norm(v1)
        v2 = -beta[j+1]*v2
        v1 = v1/beta[j+1]
    end

    #final step
    v2 = v2 + hamil*v1
    alpha[steps] = real(v1'*v2)

    return alpha, beta
end

function frac(n, z, alpha, beta)
    if n < length(alpha)
        return beta[n]^2/(z-alpha[n]-frac(n+1,z,alpha,beta))
    else
        return complex(0.0,0.0)
    end
end

function opcon(hamil,steps,ws,δ,nsite,nelec,N_c)
    eg, vg = lanczos(hamil,steps)
    J = make_Current(nsite,nelec,N_c)
    v0 = J*vg
    nm = norm(v0)
    v0 = v0/nm
    alpha, beta = lanczos(hamil,v0,steps)
    nw = length(ws)
    spec = zeros(Float64,nw)

    for i=1:nw
        z = ws[i] + δ*1im + eg
        X = frac(1,z,alpha,beta)
        spec[i] = -nm^2*imag(X)
    end

    return spec
end

function opcon(hamil,steps,ws,δ,nsite,nelec,ndown,N_c)
    eg, vg = lanczos(hamil,steps)
    J = make_Current(nsite,nelec,ndown,N_c)
    v0 = J*vg
    nm = norm(v0)
    v0 = v0/nm
    alpha, beta = lanczos(hamil,v0,steps)
    nw = length(ws)
    spec = zeros(Float64,nw)

    for i=1:nw
        z = ws[i] + δ*1im + eg
        X = frac(1,z,alpha,beta)
        spec[i] = -nm^2*imag(X)
    end

    return spec
end