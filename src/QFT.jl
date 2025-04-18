include(joinpath(@__DIR__,"General.jl"))

function QFT(L::Int64)
    d::Int64 = 2^L
    w = exp(2*pi*im/d)
    
    M = Array{ComplexF64}(undef, d, d)

    for j in 1:d
        for i in 1:d
            M[i, j] = w^( (i-1)*(j-1) )
        end
    end

    return M/sqrt(d)
end

function KrausMap(M, N) # Generates M matrices N x N
    Gin = randn(ComplexF64, M*N, N)
    Q, R = qr(Gin)
    Q = Matrix(Q)
    lambda = diagm([R[i,i]/abs(R[i,i]) for i in 1:N])

    Q = Q*lambda
    
    return [Q[(i-1)*N+1:i*N, 1:N] for i in 1:M]
end

function DilutedUnitary_QFT(r, L, k)
    d::Int64 = 2^L
    
    U = QFT(L)
    
    kraus = KrausMap(r, d)
    K = zeros(d*d, d*d)
    for i in 1:r
        K += kron(kraus[i], conj(kraus[i]))
    end

    return (1-k)*kron(U, conj(U)) + k*K
end
