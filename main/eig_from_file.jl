include(joinpath(@__DIR__,"../src/General.jl"))
include(joinpath(@__DIR__,"../src/CUE.jl"))
include(joinpath(@__DIR__,"../src/FF.jl"))
include(joinpath(@__DIR__,"../src/XXX_model.jl"))
include(joinpath(@__DIR__,"../src/Clifford.jl"))
include(joinpath(@__DIR__,"../src/QFT.jl"))

function main()

    L::Int64 = 4
    r::Int64 = 5

    d::Int64 = 2^L

    dirU = string("./../files/mat_Clif_L", L, ".dat") # substitute for desired model
    U = readdlm(dirU, ',', ComplexF64)

    dirK = string("./../files/mat_Kraus_L", L, "_r", r, ".dat")
    K = readdlm(dirK, ',', ComplexF64)

    vec_k = Array(0.01:0.01:1.)
    N = length(vec_k)

    eigvalues = Matrix{ComplexF64}(undef, 3*N, d*d)


    for i in 1:N

        k = vec_k[i]
        delta_k = k/1000

        map_1 = (1.0-k+delta_k)*U + (k-delta_k)*K
        map_2 = (1.0-k)*U + k*K
        map_3 = (1.0-k-delta_k)U + (k+delta_k)*K

        eig_1 = eigvals(map_1)
        eig_2 = eigvals(map_2)
        eig_3 = eigvals(map_3)

        eigvalues[1 + 3*(i-1), :] = eig_1
        eigvalues[2 + 3*(i-1), :] = eig_2
        eigvalues[3 + 3*(i-1), :] = eig_3

    end

    dir = string("./../files/eig_Clif_L", L, "_r", r, ".dat")
    writedlm(dir, eigvalues, ',')
end

main()
