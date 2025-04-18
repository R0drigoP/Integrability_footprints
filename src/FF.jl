include(joinpath(@__DIR__,"General.jl"))

function string_operator(state, k, l) #counts numbers of 1's
    return sum(state[k:l-1]);
end

function hilbert_base(N) #base for 2 qubits {|00>, |01>, |10>, |11>}
    base = Array{Int}(undef, 2^N, N);
    for i in 0:2^N-1
        #base[2^N-i,:] = digits(i, base=2, pad=N)|>reverse
        base[i+1,:] = digits(i, base=2, pad=N)|>reverse
    end
    
    return base
end

# 0 -> + || 1 -> -
function Hmanybody(Hnambu, N)
    base_Hilbert = hilbert_base(N)

    h = Hnambu[1:N, 1:N]
    Delta = Hnambu[1:N, N+1:end]
    d = 2^N

    base_matrix = LinearAlgebra.I(d)
    
    hmany_body = zeros(ComplexF64, d, d)
    
    for i in 1:d
        state = base_Hilbert[i,:]
        
        res = zeros(ComplexF64, d)
        diag_term = 0

        for l in 1:N #O Nadir tem isto ao contrário no código dele???
            if state[l] == 1
                diag_term += h[l, l]
            else
                diag_term -= h[l, l]
            end
            for k in 1:l-1
                p = string_operator(state, k+1, l);
                newstate = copy(state)
                if state[k] == 0
                    newstate[k] = 1
                    if state[l] == 1
                        newstate[l] = 0
                        indexket = findfirst(x -> all(x .== newstate), eachrow(base_Hilbert));
                        #res -= 2 * h[k, l] * (-1)^p * base_matrix[:, indexket]
                        res += 2 * h[k, l] * (-1)^p * base_matrix[:, indexket]
                    else
                        newstate[l] = 1
                        indexket = findfirst(x -> all(x .== newstate), eachrow(base_Hilbert));
                        res += 2 * Delta[k, l] * (-1)^p * base_matrix[:, indexket]
                    end
                else
                    newstate[k] = 0
                    if state[l] == 0
                        newstate[l] = 1
                        indexket = findfirst(x -> all(x .== newstate), eachrow(base_Hilbert));
                        res += 2 * conj(h[k, l]) * (-1)^p * base_matrix[:, indexket]
                    else
                        newstate[l] = 0
                        indexket = findfirst(x -> all(x .== newstate), eachrow(base_Hilbert));
                        #res -= 2 * conj(Delta[k, l]) * (-1)^p * base_matrix[:, indexket]
                        res += 2 * conj(Delta[k, l]) * (-1)^p * base_matrix[:, indexket]
                    end
                end
            end
        end
        
        res += diag_term * base_matrix[:, i]
        
        for j in 1:d
            hmany_body[i, j] = inner_product(base_matrix[:, j], res)
        end
    end

    return hmany_body
end

function inner_product(x::Vector, y::Vector)
    if length(x) != length(y)
        ArgumentError("Vectors must have the same length to compute inner product")
    end

    sum = 0

    for i in 1:length(x)
        sum += x[i]*y[i]
    end

    return sum
end

# Orthogonal
function rand_O(N)::Matrix{Float64}
    A = randn(Float64, N, N)
    
    Q, R = qr(A)

    lambda = diagm([R[i,i]/abs(R[i,i]) for i in 1:N])

    return Q*lambda
end

# Special Orthogonal
function rand_SO(N)::Matrix{Float64}
    U = rand_O(N)

    if det(U) < 0
        U[:, 1] = -U[:, 1]
    end
    return U
end

function rand_COE(N)::Matrix{ComplexF64}
    U = rand_CUE(N)
    
    return transpose(U)*U
end


function V_matrix(n)::Matrix{ComplexF64}
    V = BlockArray{ComplexF64}(undef_blocks, [n, n], [n, n])

    setblock!(V, Matrix(I, n, n), 1, 1)
    setblock!(V, Matrix(I, n, n), 1, 2)
    setblock!(V, -im*Matrix(I, n, n), 2, 1)
    setblock!(V, im*Matrix(I, n, n), 2, 2)
    
    
    return 1/sqrt(2)*V
    
end

function Hsinglebody(N)

    H = Array{ComplexF64}(undef, 2*N, 2*N)
    #V = V_matrix(N)
    
    u = rand_SO(2*N)
    H = im*log(u)
    
    return H

end
    

# N is the number of fermions. The matrix size is 2^N
function rand_HFF(N)

    V = V_matrix(N)
    Hnambu = V'*Hsinglebody(N)*V

    return Hmanybody(Hnambu, N)
end

function rand_UFF(N)
    H_FF = rand_HFF(N)

    #return exp(-im/2 * H_FF)
    return exp(-im * H_FF)
end

function DilutedUnitaryFF(M, L, p)
    N::Int64 = 2^L
    
    kraus = KrausMap(M, N)
    U = rand_UFF(L)

    K = zeros(N*N, N*N)
    for i in 1:M
        K += kron(kraus[i], conj(kraus[i]))
    end
    
    return (1-p)*kron(U, conj(U)) + p*K
end
