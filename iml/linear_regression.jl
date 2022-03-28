using LinearAlgebra
using CairoMakie

function SGD(x, y)
    a = 0.0
    b = 0.0

    l = Inf64

    n = length(x)

    e = 0
    i = 1

    avg_l = 0
    p_avg_l = 0

    while true
        L(a, b) = (a .* x[i] .+ b - y[i])^2

        ga = 2 * (a .* x[i] .+ b - y[i]) * (x[i])
        gb = 2 * (a .* x[i] .+ b - y[i])

        a -= 0.1 * exp(-0.1 * e) * ga
        b -= 0.1 * exp(-0.1 * e) * gb

        l = ga^2 + gb^2
        avg_l += l

        if i == n
            e += 1

            avg_l /= n

            println("y ~ $a * x + $b : ΔE[∇l] = $((avg_l-p_avg_l))")

            if e > 10 && (p_avg_l - avg_l) < 1e-12
                println("epochs = $e, E[∇l] = $avg_l")
                break
            end

            p_avg_l = avg_l
            avg_l = 0
        end

        i = mod1(i + 1, n)
    end



    (a, b)
end

