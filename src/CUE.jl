include(joinpath(@__DIR__,"General.jl"))

function DilutedUnitary(M, N, p=0, alpha=0)
    U = rand_CUE(N)
    
    if p == 0 
        if alpha == 0
            return kron(U, conj(U))
        end
            
        H = rand_GUE(N)
        noise = exp(-im*alpha*H)   
        
        return kron(noise*U, conj(noise*U)) 
    end

    kraus = KrausMap(M, N)
    K = zeros(N*N, N*N)
    for i in 1:M
        K += kron(kraus[i], conj(kraus[i]))
    end

    if alpha == 0
        return (1-p)*kron(U, conj(U)) + p*K
    end
        
    H = rand_GUE(N)
    noise = exp(-im*alpha*H)   
        
    return (1-p)*kron(noise*U, conj(noise*U)) + p*K
end