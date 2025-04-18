include(joinpath(@__DIR__,"../src/General.jl"))
include(joinpath(@__DIR__,"../src/CUE.jl"))
include(joinpath(@__DIR__,"../src/FF.jl"))
include(joinpath(@__DIR__,"../src/XXX_model.jl"))
include(joinpath(@__DIR__,"../src/Clifford.jl"))
include(joinpath(@__DIR__,"../src/QFT.jl"))

function main()

    L::Int64 = 4
    r::Vector{Int64} = [5, 2^(2*L-1), 4^L - 1] # low rank, half rank, max rank

    d::Int64 = 2^L

    U = rand_CUE(d)
    U1 = kron(U, conj(U))

    U = rand_Clifford(L)
    U2 = kron(U, conj(U))

    U = rand_UFF(L)
    U3 = kron(U, conj(U))

    U = QFT(L)
    U4 = kron(U, conj(U))

    delta = tan(rand(Uniform(0, pi)))
    U = Unitary_BW_XXXmodel(L, delta, L)
    U5 = kron(U, conj(U))
    
    kraus1 = KrausMap(r[1], d)
    K1 = zeros(d*d, d*d)
    for i in 1:r[1]
        K1 += kron(kraus1[i], conj(kraus1[i]))
    end

    kraus2 = KrausMap(r[2], d)
    K2 = zeros(d*d, d*d)
    for i in 1:r[2]
        K2 += kron(kraus2[i], conj(kraus2[i]))
    end

    kraus3 = KrausMap(r[3], d)
    K3 = zeros(d*d, d*d)
    for i in 1:r[3]
        K3 += kron(kraus3[i], conj(kraus3[i]))
    end
  
    dir = string("./../files/mat_CUE_L", L, ".dat")
    writedlm(dir, U1, ',')

    dir = string("./../files/mat_Clif_L", L, ".dat")
    writedlm(dir, U2, ',')

    dir = string("./../files/mat_FF_L", L, ".dat")
    writedlm(dir, U3, ',')

    dir = string("./../files/mat_QFT_L", L, ".dat")
    writedlm(dir, U4, ',')
 
    dir = string("./../files/mat_X_L", L, ".dat")
    writedlm(dir, U5, ',')

    dir = string("./../files/mat_Kraus_L", L, "_r", r[1], ".dat")
    writedlm(dir, K1, ',')

    dir = string("./../files/mat_Kraus_L", L, "_r", r[2], ".dat")
    writedlm(dir, K2, ',')

    dir = string("./../files/mat_Kraus_L", L, "_r", r[3], ".dat")
    writedlm(dir, K3, ',')
    
end

main()
