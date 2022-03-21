using LinearAlgebra
using GLMakie
using KrylovKit

s_x = [0 1; 
       1 0]
s_y = [0 1im; 
      -1im 0]
s_z = [1 0;
       0 -1]

eye_(N) = Matrix{Float64}(I, N, N)

function sigma_zzi(N, i)
    if i == N
        kron(s_z, eye_(2^(N - 2)), s_z)
    else
        kron(eye_(2^(N - i - 1)), s_z, s_z, eye_(2^(i - 1)))
    end
end

function sigma_xi(N, i)
    kron(eye_(2^(N - i)), s_x, eye_(2^(i - 1)))
end 

function dense_ising(N, h_x)
    sum(sigma_zzi(N, i) for i in 1:N) + h_x * sum(sigma_xi(N, i) for i in 1:N)
end

function taskA()
    N = 10
    N_h = 40
    N_E = 10

    hs = LinRange(-2, 2, N_h)
    Es = zeros(N_h, N_E)

    for i in 1:N_h
        H = Hermitian(dense_ising(N, hs[i]))

        Es[i, :] .= eigvals(H, 1:N_E + 2)[begin:N_E]
        Es[i, :] .-= Es[i, 1]
    end

    println("Plotting...")

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="h_x", ylabel="E")

    for i in 1:N_E
        lines!(hs, Es[:, i], label="E_$i")
    end

    axislegend()
    fig
end

function testDense(N)
    H = Hermitian(dense_ising(N, 1.0))
    eigvals(H)
end

## PART B

function H_z!(N, u, v)
    for s in 0:2^N - 1
        v[s + 1] = u[s + 1] * (2 * count_ones(~(s ⊻ (s >> 1 | (s & 1) << (N - 1))) & (2^N - 1)) - N)
    end
end

function H_x!(N, h_x, u, v)
    for s in 0:2^N - 1
        for i in 0:N - 1
            v[(s ⊻ (1 << i)) + 1] += u[s + 1] * h_x
        end
    end
end

function H!(N, h_x, u, v)
    H_z!(N, u, v)
    H_x!(N, h_x, u, v)
end

function testH(N=5)
    u = zeros(2^N)
    H = zeros(2^N, 2^N)

    h_x = 0.5

    for i in 1:2^N
        u[i] = 1
        H!(N, h_x, u, @view H[:, i])
        u[i] = 0
    end

    H - dense_ising(N, h_x)
end

function taskB()

    N = 4
    N_h = 100
    N_E = 16
    
    hs = LinRange(-2, 2, N_h)
    Es = zeros(N_h, N_E)

    for i in 1:N_h

        
        function H_op(u) 
            v = zeros(2^N) 
            H!(N, hs[i], u, v)
            v
        end

        eigs, vecs, info = eigsolve(H_op, rand(2^N), N_E, :SR, Lanczos())

        Es[i, :] .= eigs[begin:N_E]
    end

    println("Plotting...")

    fig = Figure()
    ax = Axis(fig[1,1], xlabel="h_x", ylabel="E")

    for i in 1:N_E
        lines!(hs, Es[:, i], label="E_$i")
    end

    axislegend()
    fig
end

function test(N)
    function H_op(u) 
        v = zeros(2^N) 
        H!(N, 1.0, u, v)
        v
    end

    eigs, vecs, info = eigsolve(H_op, rand(2^N), 16, :SR, Lanczos())

    eigs
end

