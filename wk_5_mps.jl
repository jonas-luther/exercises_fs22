

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

function main()

    N = 4
    T = zeros((ones(Int64, N).*2)...)

    T[ones(Int64, N)...] = 1
    T[(ones(Int64, N).*2)...] = 1

    T_r = reshape(T, 2, 2^(N-1))

    u,s,v = svd(T_r)
    A1 = reshape(u, 2, 1, 2)
    display(A1[1,:,:])
    display(A1[2,:,:])

    phi = diagm(s)*v'

    T_r = reshape(phi, 4, 4)
    u,s,v = svd(T_r)
    A2 = reshape(u, 2, 2, 4)

    display(A2[1,:,:])
    display(A2[2,:,:])

    phi = diagm(s)*v'

    T_r = reshape(phi, 8, 2)
    u,s,v = svd(T_r)
    A3 = reshape(u, 2, 4, 2)

    display(A3[1,:,:])
    display(A3[2,:,:])

    phi = diagm(s)*v'

    T_r = reshape(phi, 4, 1)

    u,s,v = svd(T_r)

    A4 = reshape(u, 2, 2, 1)

    display(A4[1,:,:])
    display(A4[2,:,:])

 

    display(reshape(T, 4, 4))
    display(reshape(T_rec, 4, 4))
    


   
end

function MPS(T)
    N = ndims(T)
    d = size(T)[1]

    train = []

    chi = 1

    for i in 1:N
        
        T_r = reshape(T, d * chi, d^(N-i))

        u,s,v = svd(T_r)
        chi1 = rank(T_r)

        A = reshape(u[:,1:chi1], d, chi, chi1)
        chi = chi1

        println("A_$i σ=1: $(A[1,:,:])")
        println("A_$i σ=2: $(A[2,:,:])")

        T = Diagonal(s[1:chi1])*v[:,1:chi1]'

        push!(train, A)
    end

    train[end] .*= T

    train
end

function eval_MPS(train, indices)

    A = train[1][indices[1], :, :]

    for i in 2:length(indices)
        A = A * train[i][indices[i], :, :]
    end

    A
end

function MPS_to_T(train, d, N)
    T_rec = zeros(fill(d, N)...)
    for ind in CartesianIndices(T_rec)
        println(ind)
        T_rec[ind] = eval_MPS(train, ind)[1]
    end
    T_rec
end