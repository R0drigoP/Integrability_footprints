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

# Single qb gates acting on qubit t with total of L qubits
function H(t, L)
    Had = 1/sqrt(2) * [1 1; 1 -1]
    
    return kron(Matrix(I, 2^(t-1), 2^(t-1)), Had, Matrix(I, 2^(L-t), 2^(L-t)))
end

function P(t, L)
    Phase = [1 0; 0 im]
    
    return kron(Matrix(I, 2^(t-1), 2^(t-1)), Phase, Matrix(I, 2^(L-t), 2^(L-t)))
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
#=========================================#

function empty_tableau(L)
    return Array{Int8}(undef, 2, 2*L)
end

# random sample of tableau = tab from qubit n to qubit L
function sample_rows(tab, n, L)
    for i in 1:2
        for j in n:L
            tab[i, j] = rand([0, 1])
            tab[i, L+j] = rand([0, 1])
        end
    end
end # return void, want to alter tab directly

function check_anti_commute_tab(tab, n, L)
    commute = 0

    for i in n:L
        x1 = tab[1, i]
        z1 = tab[1, L+i]

        x2 = tab[2, i]
        z2 = tab[2, L+i]

        if (x1 == z1 == 0) || (x2 == z2 == 0) || (x1 == x2 && z1 == z2) 
            commute += 1
        end
    end
        
    return (L-n+1-commute)%2 == 1
end

function sample_anti_commute_rows(tab, n, L)
    while true
        sample_rows(tab, n, L)
        check_anti_commute_tab(tab, n, L) && break
    end
end

function sweep(tab, n, L)
    M = Matrix(I, 2^L, 2^L)

    # step 1
    for i in n:L
        if tab[1, L+i] == 1
            if tab[1, i] == 0
                M = H(i, L)*M
                swapcols!(tab, i, L+i)
            else
                M = P(i, L)*M
                tab[1, L+i] = 0
                tab[2, L+i] = (tab[2, i] + tab[2, L+i])%2
            end
        end
    end


    # step 2
    J = [i for i in n:L if tab[1, i] == 1]
    N_J = length(J)
    for i in 1:trunc(Int, N_J/2) + N_J%2
        used = []
        for i in 1:2:length(J) - 1
            M = CNOT(J[i], J[i+1], L)*M
            tab[1, J[i+1]] = (tab[1, J[i]] + tab[1, J[i+1]] )%2
            tab[2, J[i+1]] = (tab[2, J[i]] + tab[2, J[i+1]] )%2
            
            tab[1, L+J[i]] = (tab[1, L+J[i]] + tab[1, L+J[i+1]] )%2
            tab[2, L+J[i]] = (tab[2, L+J[i]] + tab[2, L+J[i+1]] )%2
            
            append!(used, i+1)
        end
        deleteat!(J, used)
    end
    
    # step 3
    if J != [n]
        M = SWAP(n, J[1], L)*M

        # CNOT(n, J[1])
        tab[1, J[1]] = (tab[1, n] + tab[1, J[1]] )%2
        tab[2, J[1]] = (tab[2, n] + tab[2, J[1]] )%2
        
        tab[1, L+n] = (tab[1, L+n] + tab[1, L+J[1]] )%2
        tab[2, L+n] = (tab[2, L+n] + tab[2, L+J[1]] )%2

        # CNOT(J[1], n)
        tab[1, n] = (tab[1, n] + tab[1, J[1]] )%2
        tab[2, n] = (tab[2, n] + tab[2, J[1]] )%2
        
        tab[1, L+J[1]] = (tab[1, L+n] + tab[1, L+J[1]] )%2
        tab[2, L+J[1]] = (tab[2, L+n] + tab[2, L+J[1]] )%2

        # CNOT(n, J[1])
        tab[1, J[1]] = (tab[1, n] + tab[1, J[1]] )%2
        tab[2, J[1]] = (tab[2, n] + tab[2, J[1]] )%2
        
        tab[1, L+n] = (tab[1, L+n] + tab[1, L+J[1]] )%2
        tab[2, L+n] = (tab[2, L+n] + tab[2, L+J[1]] )%2
    end

    
    #step 4
    if !(tab[2, n] == 0 && tab[2, L+n] == 1 && tab[2, n+1:L] == zeros(L-n) && tab[2, L+n+1:2*L] == zeros(L-n))
        M = H(n, L)*M
        swapcols!(tab, n, L+n)

        #Repeat step 1 for second row
        for i in n:L
            if tab[2, L+i] == 1
                if tab[2, i] == 0
                    M = H(i, L)*M
                    swapcols!(tab, i, L+i)
                else
                    M = P(i, L)*M
                    tab[1, L+i] = (tab[1, i] + tab[1, L+i])%2
                    tab[2, L+i] = 0
                end
            end
        end

        # Repeat step 2 for second row
        J = [i for i in n:L if tab[2, i] == 1]
        N_J = length(J)
        for i in 1:trunc(Int, N_J/2) + N_J%2
            used = []
            for i in 1:2:length(J) - 1
                M = CNOT(J[i], J[i+1], L)*M
                tab[1, J[i+1]] = (tab[1, J[i]] + tab[1, J[i+1]] )%2
                tab[2, J[i+1]] = (tab[2, J[i]] + tab[2, J[i+1]] )%2
                
                tab[1, L+J[i]] = (tab[1, L+J[i]] + tab[1, L+J[i+1]] )%2
                tab[2, L+J[i]] = (tab[2, L+J[i]] + tab[2, L+J[i+1]] )%2
                
                append!(used, i+1)
            end
            deleteat!(J, used)
        end
    
         # Repeat step 3 for second row
        if J != [n]
            M = SWAP(n, J[1], L)*M
    
            # CNOT(n, J[1])
            tab[1, J[1]] = (tab[1, n] + tab[1, J[1]] )%2
            tab[2, J[1]] = (tab[2, n] + tab[2, J[1]] )%2
            
            tab[1, L+n] = (tab[1, L+n] + tab[1, L+J[1]] )%2
            tab[2, L+n] = (tab[2, L+n] + tab[2, L+J[1]] )%2
    
            # CNOT(J[1], n)
            tab[1, n] = (tab[1, n] + tab[1, J[1]] )%2
            tab[2, n] = (tab[2, n] + tab[2, J[1]] )%2
            
            tab[1, L+J[1]] = (tab[1, L+n] + tab[1, L+J[1]] )%2
            tab[2, L+J[1]] = (tab[2, L+n] + tab[2, L+J[1]] )%2
    
            # CNOT(n, J[1])
            tab[1, J[1]] = (tab[1, n] + tab[1, J[1]] )%2
            tab[2, J[1]] = (tab[2, n] + tab[2, J[1]] )%2
            
            tab[1, L+n] = (tab[1, L+n] + tab[1, L+J[1]] )%2
            tab[2, L+n] = (tab[2, L+n] + tab[2, L+J[1]] )%2
        end
        
        #Repeat step 4 for second row
        if !(tab[2, n] == 0 && tab[2, L+n] == 1 && tab[2, n+1:L] == zeros(L-n) && tab[2, L+n+1:2*L] == zeros(L-n))
            M = H(n, L)*M
            swapcols!(tab, n, L+n)
        end
    end
    
    #step 5
    sa = rand([0, 1])
    sb = rand([0, 1])

    if sa == 0 && sb == 1
        M = X(n, L)*M
    elseif sa == 1 
        if sb == 1
            M = Y(n, L)*M
        else
            M = Z(n, L)*M
        end
    end
    
    return M 
end

# L = number of qubits
function rand_Clifford(L)
    M = I

    tab = empty_tableau(L)

    for i in 1:L
        sample_anti_commute_rows(tab, i, L)
        M *= sweep(tab, i, L)
    end

    return M
end

function swapcols!(X::AbstractMatrix, i::Integer, j::Integer)
    @inbounds for k = 1:size(X,1)
        X[k,i], X[k,j] = X[k,j], X[k,i]
    end
end

function DilutedUnitaryClifford(M, L, p=0, alpha=0)
    N::Int64 = 2^L
    U = rand_Clifford(L)
    
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