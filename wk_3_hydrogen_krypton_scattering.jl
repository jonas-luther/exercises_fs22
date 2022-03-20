using OffsetArrays
using GLMakie
using SpecialFunctions

function numerov_step(psi, psi_prev, x, dx, K)
    x_prev = x - dx
    x_next = x + dx

    c = dx^2 / 12

    a = 2(1 - 5 * c * K(x)) * psi - (1 + c * K(x_prev)) * psi_prev

    a / (1 + c * K(x_next))
end

function numerov_integrate(N, x0, x1, K, psi0, normalize=true)

    dx = (x1 - x0) / N

    psi = OffsetArray{ComplexF64}(zeros(N + 1), 0:N)
    psi[0] = psi0(x0)
    psi[1] = psi0(x0 + dx)

    for i in 1:(N - 1)
        psi[i + 1] = numerov_step(psi[i], psi[i - 1], x0 + i * dx, dx, K)
    end

    if normalize
        psi_norm =  sqrt(sum(psi[i]^2 * dx for i in 0:N))
        psi ./= psi_norm
    end

    LinRange(x0, x1, N + 1), psi[0:N]
end


bj(ν, x) = √(π / 2x) * besselj(ν + 1 / 2, x)
bn(ν, x) = √(π / 2x) * bessely(ν + 1 / 2, x)


function main()

    N = 1000

    A = 6.12
    epsilon = 5.9
    C = sqrt(A * epsilon) / 5

    r_min = 0.5
    r_max = 5

    function delta(E, l)

        K(r) = (E - epsilon * (r^-12 - 2 * r^-6)) * A - l * (l + 1) * r^-2
        psi0(r) = exp(-C * r^-5)
    
        rs, u_l = numerov_integrate(N, r_min, r_max, K, psi0, true)

        k0 = sqrt(A * E)

        r1 = rs[N]
        r2 = rs[N + 1]

        u1 = u_l[N]
        u2 = u_l[N + 1]

        K_ =  r1 * u2 / (r2 * u1)

        atan((K_ * bj(l, k0 * r1) - bj(l, k0 * r2)) / (K_ * bn(l, k0 * r1) - bn(l, k0 * r2)))
    end

    function total_scatter(E, l_max)
        4 * pi / (A * E) * sum([(2l + 1) * sin(delta(E, l))^2 for l in 0:l_max])
    end

    Es = 0.1:0.02:3.5
    lines(Es, [real(total_scatter(E, 10)) for E in Es], axis=(xticks = LinearTicks(8), xlabel = "E [meV]", ylabel = "σ/σ²"), label="l_max = 10")

    axislegend()
end


main()
current_figure()
