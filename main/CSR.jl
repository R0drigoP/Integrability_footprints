include(joinpath(@__DIR__,"../src/General.jl"))
include(joinpath(@__DIR__,"../src/CUE.jl"))
include(joinpath(@__DIR__,"../src/FF.jl"))
include(joinpath(@__DIR__,"../src/XXX_model.jl"))
include(joinpath(@__DIR__,"../src/Clifford.jl"))
include(joinpath(@__DIR__,"../src/QFT.jl"))

function main()

    L::Int64 = 6
    It::Int64 = 50

    d::Int64 = 2^L
    r::Int64 = 5

    k_min::Float64 = 0.02
    k_step::Float64 = 0.02
    k_max::Float64 = 1.

    vec_k::Array{Float64} = Array(k_min:k_step:k_max)
    N::Int64 = length(vec_k)
    N_eig::Int64 = d*d-1

    CSR::Array{ComplexF64} = Array{Float64}(undef, It*N_eig, N)
    av_radius::Array{Float64} = zeros(N)

    for iter in 1:It
	#Choose ensemble	

        #U = rand_CUE(d)
        #U = rand_Clifford(L)
        U = rand_UFF(L)
        #delta = tan(rand(Uniform(0, pi)))
        #U = Unitary_BW_XXXmodel(L, delta, L)

        kraus = KrausMap(r, d)

        for i in 1:N

            k = vec_k[i]
            map = DilutedUnitaryFromMatrices(U, kraus, k)
            eig = eigvals(map)
            pop!(eig)

	    CSR[1+(iter-1)*N_eig:iter*N_eig, i] = complex_spacing_ratio(eig)

        end
    end
	
    dir = string("./../files/FF_CSR_L", L, "_r", r, "_It", It, ".dat")
    writedlm(dir, CSR, ',')

end

main()