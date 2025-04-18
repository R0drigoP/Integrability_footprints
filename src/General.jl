using LinearAlgebra
using Statistics
using Random
using Kronecker
using BlockArrays
using DelimitedFiles
#using PoissonRandom
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

function plot_cumulative_distribution(data::Vector{<:Real}; bins::Int=50)
    # Create a histogram
    h = fit(Histogram, data, nbins=bins)

    # Compute the cumulative distribution
    cumulative_counts = cumsum(h.weights)
    cumulative_distribution = cumulative_counts ./ sum(h.weights)

    # Compute bin centers
    bin_centers = h.edges[1][1:end-1] .+ diff(h.edges[1]) ./ 2

    return bin_centers, cumulative_distribution
end

#k = 1 -> nearest neighbor
function kth_NN(data, k=1)
    N = length(data)
    k_NN::Vector{Int64} = zeros(N)
    
    points = hcat(real(data), imag(data))'
    
    kdtree = KDTree(points)
    idxs, dists = knn(kdtree, points, k+1, true)

    #display(idxs[1][1])
    
    for i in 1:N
        k_NN[i] = idxs[i][k+1]
    end
    
    return k_NN
end

# Nearest neighbour of data1 in data2
function find_nearest_neighbors_complex(data1::Vector{<:Complex}, data2::Vector{<:Complex})
    # Convert complex numbers to 2D points (real and imaginary parts as Float64)
    points1 = hcat(real(data1), imag(data1))' .|> Float64   # Each column is a 2D point
    points2 = hcat(real(data2), imag(data2))' .|> Float64

    # Build a KDTree from points2
    tree = KDTree(points2)
    
    # Prepare arrays to store nearest neighbor indices and distances
    nearest_indices = Array{Int64, 1}(undef, length(data1))
    nearest_distances = Array{Float64, 1}(undef, length(data1))

    # For each point in data1, find the nearest neighbor in data2
    for i in 1:length(data1)
        query_point = points1[:, i]
        idx, dist = knn(tree, query_point, 1)  # find the nearest neighbor
        nearest_indices[i] = idx[1]            # store the index of the nearest neighbor
        nearest_distances[i] = dist[1]         # store the distance to the nearest neighbor
    end

    return nearest_indices, nearest_distances
end

# Kolmogoriv-Smirov test. hyp -> hypothesis, x, y -> coordinates of the test function
function KS_test(hyp::Function, x, y)
    return maximum([abs(hyp(x[i]) - y[i]) for i in 1:length(x)])
end

function exp_pdf(x, N=1, a=0, b=1)
    lambda = N/(b-a)

    return lambda*exp(-lambda*x)
end

# cumulative
function exp_cdf(x, N=1, a=0, b=1)
    lambda = N/(b-a)

    return 1-exp(-lambda*x)
end

function wigner_pdf(x)
    return 0.5*pi*x*exp(-0.25*pi*x*x)
end

function wigner_cdf(x)
    return 1-exp(-0.25*pi*x*x)
end

# n is order of polynomial aprox, N is number of data points
function ginibre_cdf(x, N=1000)
    res = 1

    for n in 1:N-1
        res *= exp(Taylor1(Float64, n))(x*x) * exp(-x*x)
    end

    return 1-res
end

# z0 is target; zj is vector with all data; N is sample size
function rho_unfolded(z0, zj, N, sigma)
    rho::Float64 = 0

    for i in 1:N
        #x = exp( -abs(z0 - zj[i])^2/(2*sigma*sigma) )
        rho += exp( -abs(z0 - zj[i])^2/(2*sigma*sigma) )
        #rho += x
        #display(x)
    end
    
    #rho::Float64 = sum(exp.(-1/(2*sigma*sigma) * abs.(z0 .- zj).^2))

    return rho/(2*pi*sigma*sigma*N)
    #return rho/(sqrt(2*pi)*sigma)
end

function S_unfolded(S_rho, N; nbins=100)
    hist = fit(Histogram, S_rho, nbins=nbins)

    S = midpoints(hist.edges[1])
    p_S = hist.weights
    
    norm_p_S = p_S / sum(p_S) #normalized bins

    first_moment = sum(S .* norm_p_S)

    display(first_moment)
    
    norm_pS = norm_p_S * first_moment #normalized with first moment equal to 1
    norm_S = S*first_moment

    return norm_S, norm_pS
end

function heaviside(x)
   0.5 * (sign(x) + 1)
end

function rho_ginibre(z, M)
    return 1/(M*pi)* abs(z)^(2*(1-M)/M)*heaviside(1-abs(z))
end

function p_ginibre(s, N=10)
    res1 = 1
    res2 = 0

    S = 1.14*s
    
    for k in 1:N
        res1 *= gamma(1+k, S*S)/gamma(1+k)
        res2 += 2*S^(2*k+1)*exp(-S*S)/gamma(1+k, S*S)
    end

    return res1*res2*1.14
end

function rand_GinUE(N)
    return randn(ComplexF64, N, N)/sqrt(N)
end

function get_histogram_points(data; nbins=100)

    # Create a histogram
    hist = fit(Histogram, data, nbins=100)
    
    # Calculate the bin centers (x values)
    bin_centers = midpoints(hist.edges[1])
    
    # Get the bin heights (y values)
    bin_heights = hist.weights

    # Calculate the bin widths
    bin_widths = diff(hist.edges[1])
    
    # Normalize the bin heights to represent a PDF
    pdf_heights = bin_heights ./ (sum(bin_heights .* bin_widths))
    
    return bin_centers, pdf_heights
end

# Function to compute the first moment
function first_moment(x_data, rho_data)
    # Interpolate the PDF
    itp = extrapolate(interpolate((x_data,), rho_data, Gridded(Linear())), Flat())

    # Integrate x * ρ(x) over the range of the PDF to compute the first moment
    M_1, _ = quadgk(x -> x * itp(x), minimum(x_data), maximum(x_data))
    
    return M_1
end

# Function to transform the PDF as ρ(x) -> M_1 * ρ(M_1 * x)
function transform_pdf(x_data, rho_data)
    M_1 = first_moment(x_data, rho_data)
    
    # Interpolate the original PDF
    itp = extrapolate(interpolate((x_data,), rho_data, Gridded(Linear())), Flat())

    # Define the transformed PDF
    function transformed_density(x)
        M_1 * itp(M_1 * x)
    end

    # Generate a new range of x values after transformation
    x_transformed = range(0, stop=3, length=500)
    rho_transformed = transformed_density.(x_transformed)
    
    return x_transformed, rho_transformed
end
