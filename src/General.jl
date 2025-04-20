using LinearAlgebra
using Statistics
using Random
using Kronecker
using BlockArrays
using DelimitedFiles
using RandomNumbers
using Distributions
using StatsBase
using NearestNeighbors
using TaylorSeries
using SpecialFunctions
using NumericalIntegration
using QuadGK
#using Interpolations

function rand_CUE(N)::Matrix{ComplexF64}
    A = randn(ComplexF64, N, N)
    
    Q, R = qr(A)
    lambda = diagm([R[i,i]/abs(R[i,i]) for i in 1:N])
    
    return Q*lambda
end

function rand_GUE(n)::Matrix{ComplexF64}
    M = Array{ComplexF64}(undef, n, n)
    for i in 1:n
        M[i, i] = randn(Float64)
        for j in i+1:n
            M[i, j] = randn(ComplexF64)
            M[j, i] = conj(M[i, j])
        end
    end
    
    return M
end

function KrausMap(M, N) # Generates M matrices N x N
    Gin = randn(ComplexF64, M*N, N)
    Q, R = qr(Gin)
    Q = Matrix(Q)

    return [Q[(i-1)*N+1:i*N, 1:N] for i in 1:M]
end

function DilutedUnitaryFromMatrices(U, kraus, k)
    N::Int64 = length(U[1,:])
    rank::Int64 = length(kraus)

    K = zeros(N*N, N*N)
    for i in 1:rank
        K += kron(kraus[i], conj(kraus[i]))
    end
    
    return (1-k)*kron(U, conj(U)) + k*K
end

#k = 1 -> nearest neighbor
function kth_NN(data, k=1)
    N = length(data)
    k_NN::Vector{Int64} = zeros(N)
    
    points = hcat(real(data), imag(data))'
    
    kdtree = KDTree(points)
    idxs, dists = knn(kdtree, points, k+1, true)

    
    for i in 1:N
        k_NN[i] = idxs[i][k+1]
    end
    
    return k_NN
end

function rand_GinUE(N)
    return randn(ComplexF64, N, N)/sqrt(N)
end

