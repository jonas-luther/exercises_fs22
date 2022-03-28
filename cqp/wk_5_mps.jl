

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

function MPS(T)
    N = ndims(T)
    d = size(T)[1]

    train = Vector{Array{Float64,3}}()

    chi = 1
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


function MPS_to_T(train, d, N)
    T_rec = zeros(fill(d, N)...)

    for indices in CartesianIndices(T_rec)
        A = train[1][:, indices[1], :]

        for i in 2:length(indices)
            A = A * train[i][:, indices[i], :]
        end

        T_rec[indices] = A[1]
    end

    T_rec
end

function print_MPS(train)
    for i in 1:length(train)
        println("A_$i σ=1: ")
        display(train[i][:, 1, :])
        println("A_$i σ=2: ")
        display(train[i][:, 2, :])
    end
end

function MPS_to_vidal(train, d, N)


end