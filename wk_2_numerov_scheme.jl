module Week2

using GLMakie
using OffsetArrays


function numerov(psi, psi_prev, x, dx, K)
    x_prev = x - dx
    x_next = x + dx

    c = dx^2 / 12

    a = 2(1 - 5 * c * K(x)) * psi - (1 + c * K(x_prev)) * psi_prev

    a / (1 + c * K(x_next))
end

function main()

    ħ = 1

    N = 10000

    dx = 1 / N 

    c = 10000

    function V(x)
        c * (x^2 - x)
    end

    function get_psi(E)

        m = 1

        b = 0.5 + 0.5 * sqrt(1 + 4E / c)

        b_index = Int64(floor(b ÷ dx))

        function K(x)
            2m * (E - V(x)) / ħ^2
        end

        psiL = OffsetArray{ComplexF64}(zeros(N + 3), -1:N + 1)

        psiL[0] = 1
        psiL[-1] = exp(-dx * sqrt(Complex(2 * m * E)) / ħ)

        for i in 0:N
            psiL[i + 1] = numerov(psiL[i], psiL[i - 1], i * dx, dx, K)
        end

        psiR = OffsetArray{ComplexF64}(zeros(N + 3), -1:N + 1)

        psiR[N] = 1
        psiR[N + 1] = exp(-dx * sqrt(Complex(2 * m * E)) / ħ)

        for i in N:-1:0
            psiR[i - 1] = numerov(psiR[i], psiR[i + 1], i * dx, -dx, K)
        end

        nrm = psiR[b_index] / psiL[b_index]

        derL = (psiL[b_index + 1] - psiL[b_index - 1]) / 2dx
        derR = (psiR[b_index + 1] - psiR[b_index - 1]) / 2dx

        [psiL[0:b_index - 1]; psiR[b_index:N] / nrm], derL / psiL[b_index] - derR / psiR[b_index]
    end

    fig = Figure()

    display(fig)

    Es = -0.25 * c:0.01:0.0

    ax = Axis(fig[1,1])
    s =  Slider(fig[2,1][1,1], range=Es)
    root_find_button = Button(fig[2,1][1,2])

    psi_obs = lift(s.value) do E
        psi, err = get_psi(E)

        psi[1:N + 1] / sqrt((psi[1:N + 1]' * psi[1:N + 1])) / dx
    end

    lines!(fig[1,1], 0:dx:1, lift(psi_obs) do psi
        real.(psi)
    end, label="Re Ψ")

    lines!(fig[1,1], 0:dx:1, lift(psi_obs) do psi
        imag.(psi)
    end, label="Im Ψ")

    function error(E)
        psi, err = get_psi(E)
        err
    end

    axislegend()
    
    errs = error.(Es)

    lines(fig[1,2], Es, real.(errs), label="Re err")
    lines!(fig[1,2], Es, imag.(errs), label="Im err")

    on(root_find_button.clicks) do 
        s.value[] = -50.0
    end

    ylims!(-40, 40.0)

    axislegend()

    scatter!(fig[1,2], lift(s.value) do E 
        [E]
    end, [0])

    # lines!(fig[1,1], 0:dx:1, V.(0:dx:1), label="V")

    

    fig
end
end