

module Ex5_1
function main()
    vup = 1 / sqrt(2) * [0 1]
    vdown = 1 / sqrt(2) * [1 0]
    Aup = [0 1; 0 0]
    Adown = [0 0; 1 0]
    uup = [1, 0]
    udown = [0, 1]

    vup * Adown * Aup * Adown * Aup * udown
end
end

using LinearAlgebra
using Einsum

function MPS(T)
    N = ndims(T) # dimensionality of tensor
    d = size(T)[1] # degrees of freedom per site (2 for spin-1/2)

    train = Vector{Array{Float64,3}}() # will store all A matrices

    chi = 1 # first bond dimension is always 1 (A_1 is a 1xchi2 matrix)

    for i in 1:N
        T_r = reshape(T, d * chi, d^(N - i))
        u, s, v = svd(T_r)
        chi1 = rank(T_r, atol=1e-10, rtol=1e-6)
        A = reshape(u[:, 1:chi1], chi, d, chi1) # tensor network defines index sequence !!
        chi = chi1
        T = Diagonal(s[1:chi1]) * (v[:, 1:chi1]')
        push!(train, A)
    end

    train[end] .*= T
    train
end


function dims(train)
    d = size(train[1])[2]
    N = length(train)

    d, N
end

function MPS_to_T(train)

    d, N = dims(train)

    T_rec = zeros(fill(d, N)...)

    for indices in CartesianIndices(T_rec)
        T_rec[indices] = prod(train[i][:, indices[i], :] for i in 1:N)[1]
    end

    T_rec
end

function print_MPS(train)

    d, N = dims(train)

    for i in 1:N
        for l in 1:d
            println("A_$i σ=$l: ")
            display(train[i][:, l, :])
        end
    end
end

function check_canonical(train)
    d, N = dims(train)

    left_canonical = fill(false, N)
    right_canonical = fill(false, N)

    for i in 1:N
        left_canonical[i] = norm(sum(train[i][:, l, :]' * train[i][:, l, :] for l in 1:d) - I) < 1e-6
        right_canonical[i] = norm(sum(train[i][:, l, :] * train[i][:, l, :]' for l in 1:d) - I) < 1e-6
    end

    left_canonical, right_canonical
end

function MPS_to_vidal(train)
    d, N = dims(train)

    Ms = copy(train)

    Bs = Vector{Array{Float64,3}}() # will store all B tensors

    for i in N:-1:2
        da, s, db = size(Ms[i])

        M = reshape(Ms[i], da, s * db)

        U, S, V = svd(M)

        B = reshape(V', da, s, db)

        pushfirst!(Bs, B)

        for l in 1:d
            Ms[i-1][:, l, :] .= Ms[i-1][:, l, :] * U * Diagonal(S)
        end
    end

    pushfirst!(Bs, train[1]) # will only be left normalized if MPS is normalized

    Gammas = Vector{Array{Float64,3}}() # will store all Γ tensors
    Lambdas = Vector{Diagonal{Float64}}() # will store all Λ matrices

    for i in 1:N-1
        da, s, db = size(Bs[i])

        M = reshape(Bs[i], s * da, db)

        U, S, V = svd(M)

        push!(Gammas, zeros(da, s, db))
        push!(Lambdas, Diagonal(S))

        for l in 1:s
            lam_inv = i > 1 ? inv(Lambdas[i-1]) : I

            Gammas[end][:, l, :] = lam_inv * reshape(U, da, s, db)[:, l, :]
            Bs[i+1][:, l, :] = (Diagonal(S) * V') * Bs[i+1][:, l, :]
        end
    end

    push!(Gammas, Bs[N])

    for l in 1:d
        Gammas[end][:, l, :] = inv(Lambdas[end]) * Gammas[end][:, l, :]
    end

    Lambdas, Gammas

end


function evaluate_vidal(Lambdas, Gammas, indices)
    _, N = dims(Gammas)

    train = []

    for i in 1:N
        push!(train, Gammas[i][:, indices[i], :])
        if i < N
            push!(train, Lambdas[i])
        end
    end

    prod(train)[1]
end

function vidal_to_T(Lambdas, Gammas)
    d, N = dims(Gammas)

    T_rec = zeros(fill(d, N)...)

    for indices in CartesianIndices(T_rec)
        T_rec[indices] = evaluate_vidal(Lambdas, Gammas, indices)
    end

    T_rec
end

function vidal_product(LambdasBra, GammasBra, LambdasKet, GammasKet)

    d, N = dims(GammasBra)

    ketl = LambdasKet[1]
    ketg = GammasKet[1]
    bral = LambdasBra[1]
    brag = GammasBra[1]

    @einsum A[i, j] := ketl[i, i] * ketg[1, i, k] * brag[1, k, j] * bral[j, j]

    for i in 2:N-1
        ketl = LambdasKet[i]
        ketg = GammasKet[i]
        bral = LambdasBra[i]
        brag = GammasBra[i]

        
    end



    display(A)
end
