include(joinpath(@__DIR__,"General.jl"))

# CNOT gate with control c, target t, with L total qubits
function CNOT(c, t, L)
    proj_0 = [1 0; 0 0]
    proj_1 = [0 0; 0 1]
    X = [0 1; 1 0]

    if t > c
        return kron(Matrix(I, 2^(c-1), 2^(c-1)), proj_0, Matrix(I, 2^(L-c), 2^(L-c))) + 
        kron(Matrix(I, 2^(c-1), 2^(c-1)), proj_1, Matrix(I, 2^(t-c-1), 2^(t-c-1)), X, Matrix(I, 2^(L-t), 2^(L-t)))
    else
        return kron(Matrix(I, 2^(c-1), 2^(c-1)), proj_0, Matrix(I, 2^(L-c), 2^(L-c))) + 
        kron(Matrix(I, 2^(t-1), 2^(t-1)), X, Matrix(I, 2^(c-t-1), 2^(c-t-1)), proj_1, Matrix(I, 2^(L-c), 2^(L-c)))
    end
end

# SWAP gate of qubits c and t with total of L qubits
function SWAP(c, t, L)
    CNOT1 = CNOT(c, t, L)
    CNOT2 = CNOT(t, c, L)
    return CNOT1*CNOT2*CNOT1
end

function X(t, L)
    Pauli_X = [0 1; 1 0]
    
    return kron(Matrix(I, 2^(t-1), 2^(t-1)), Pauli_X, Matrix(I, 2^(L-t), 2^(L-t)))
end

function Y(t, L)
    Pauli_Y = [0 -im; im 0]
    
    return kron(Matrix(I, 2^(t-1), 2^(t-1)), Pauli_Y, Matrix(I, 2^(L-t), 2^(L-t)))
end

function Z(t, L)
    Pauli_Z = [1 0; 0 -1]
    
    return kron(Matrix(I, 2^(t-1), 2^(t-1)), Pauli_Z, Matrix(I, 2^(L-t), 2^(L-t)))
end

# Applies a given 2-qb gate U to qubits t and t+1, on a system with L qubits
function two_qb_U(U, t, L)
    return kron(Matrix(I, 2^(t-1), 2^(t-1)), U, Matrix(I, 2^(L-t-1), 2^(L-t-1)))
end

# Applies a given 2-qb gate mqp to qubits t and t+1, on a system with L qubits
function two_qb_K(K, t, L)
    return kron(Matrix(I, 4^(t-1), 4^(t-1)), K, Matrix(I, 4^(L-t-1), 4^(L-t-1)))
end

# R matrix of XXX model acting on 2 qubits
function R_XXXmodel(delta)
    P = 0.5*(I + X(1, 2)*X(2, 2) + Y(1, 2)*Y(2, 2) + Z(1, 2)*Z(2, 2) )

    return (I + im*delta*P)/(1 + delta*im)
end

# L must be even, since we are considering Periodic BC
function Unitary_BW_XXXmodel(L, delta, Tmax)
    d::Int64 = 2^L
    U = Matrix(I, d, d)
    R = R_XXXmodel(delta)

    # even part
    for q in 2:2:L-1
        U *= two_qb_U(R, q, L)
    end

    U *= SWAP(2, L, L)*two_qb_U(R, 1, L)*SWAP(2, L, L)
    
    # odd part
    for q in 1:2:L
        U *= two_qb_U(R, q, L)
    end

    return U^Tmax
end

# L is even
function kraus_BW_XXXmodel(r, L, Tmax)
    d::Int64 = 2^L

    kraus = KrausMap(r, 4)
    K = zeros(16, 16)
     for i in 1:r
        K += kron(kraus[i], conj(kraus[i]))
    end

    U = Matrix(I, d*d, d*d)

    for t in 1:Tmax    
        # odd part
        for q in 1:2:L
            U *= two_qb_K(K, q, L)
        end

        # even part
        for q in 2:2:L-1
            U *= two_qb_K(K, q, L)
        end

        U *= kron(SWAP(2, L, L), SWAP(2, L, L))*two_qb_K(K, 1, L)*kron(SWAP(2, L, L), SWAP(2, L, L))
    end

    return U 
end

function DilutedUnitaryXXXmodelIntegrable(M, L, delta, Tmax, p)
    N::Int64 = 2^L
    
    kraus = KrausMap(M, N)
    U = Unitary_BW_XXXmodel(L, delta, Tmax)

    K = zeros(N*N, N*N)
    for i in 1:M
        K += kron(kraus[i], conj(kraus[i]))
    end
    
    return (1-p)*kron(U, conj(U)) + p*K
end
