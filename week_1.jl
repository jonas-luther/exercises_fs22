
module Week1

using LinearAlgebra
using GLMakie

include("QM.jl")
using .QM

# plots magnetization in relation to inverse temperature beta
function part1()

    h_z = 1

    Ĥ = Hermitian(- h_z * s_z)
    eig = eigen(Ĥ)

    function density(beta)
        rhoZ = exp(-beta * Ĥ)

        rhoZ / tr(rhoZ)
    end

    betas = 0:0.03:3.0
    purities = [tr(density(beta)^2) for beta in betas]
    avgs_z = [real(tr(density(beta) * s_z)) for beta in betas]
    avgs_y = [real(tr(density(beta) * s_y)) for beta in betas]

    lines(betas, avgs_z, label="Mz")
    lines!(betas, avgs_y, label="My")
    lines!(betas, purities, label="purity")
    axislegend()

    current_figure()
end

# plots magnetization in relation to inverse temperature beta
function part1T()

    h_z = 1

    Ĥ = Hermitian(- h_z * s_z)
    eig = eigen(Ĥ)

    function density(beta)
        rhoZ = exp(-beta * Ĥ)

        rhoZ / tr(rhoZ)
    end

    Ts = 0.00001:0.01:100
    purities = [tr(density(1 / T)^2) for T in Ts]
    avgs_z = [real(tr(density(1 / T) * s_z)) for T in Ts]
    avgs_y = [real(tr(density(1 / T) * s_y)) for T in Ts]

    lines(Ts, avgs_z, label="Mz")
    lines!(Ts, avgs_y, label="My")
    lines!(Ts, purities, label="purity")
    axislegend()

    current_figure()
end


function two_spin_hamiltonian(lam)
    s_plus = (s_x + im * s_y) / √2
    s_minus = (s_x - im * s_y) / √2

    Hermitian(lam * (kron(s_plus, s_minus) + kron(s_minus, s_plus)) + kron(s_z, s_z))
end

function density_matrix(H, beta)
    expH = exp(-beta * H)

    expH / tr(expH)
end


# plots total magnetization of two spins in relation to inverse temperature
function part2(lam)
    H = two_spin_hamiltonian(lam)

    energies, states = eigen(H)

    density_matrix(H, 10)

    for (i, energy) in enumerate(energies)
        println("|$i> = $(states[:,i])")
        println("H|$i> = $(real(energy)) |$i>")
        println("<$i|S_z_1|$i> = $(real(states[:,i]' * s_z_1 * states[:,i]))")
        println("<$i|S_z_2|$i> = $(real(states[:,i]' * s_z_2 * states[:,i]))")
    end

    betas = 0.0:0.01:3.0

    M = [real(tr(density_matrix(H, beta) * (s_z_1 + s_z_2)^2)) for beta in betas]
    lines(betas, M)
end

# plots time evolution of psi in |+> and |-> basis (we only observe dynamics when psi(0) is not an eigenstate of H
# |+> and |-> are eigenstates of σ_y
function part3()
    plus = [1, im] / √2
    minus = [1, -im] / √2

    psi0 = plus

    psi0 /= sqrt(psi0' * psi0)

    H = 1.0 * -s_z # magnetic field in z direction

    U(t) = exp(im * t * H)

    times = 0:0.01:4π

    psis = [U(t) * psi0 for t ∈ times]

    lines(times, [abs2(plus' * psi) for psi ∈ psis], label="|<+|ψ(t)>|²", axis=(xticks = MultiplesTicks(9, π, "π"),))
    lines!(times, [abs2(minus' * psi) for psi ∈ psis], label="|<-|ψ(t)>|²")
    axislegend()

    current_figure()
end

end    