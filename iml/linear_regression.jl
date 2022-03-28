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

    Ls = Vector{Float64}()

    va = 0
    vb = 0

    while true
        L(a, b) = (a .* x[i] .+ b - y[i])^2

        ga = 2 * (a .* x[i] .+ b - y[i]) * (x[i])
        gb = 2 * (a .* x[i] .+ b - y[i])

        nu = exp(-e)

        nva = 0.5 * nu * (ga + va)
        nvb = 0.5 * nu * (gb + vb)

        a -= nva
        b -= nvb

        va = nva
        vb = nvb

        l = ga^2 + gb^2
        avg_l += l

        push!(Ls, norm(a * x .+ b - y))

        if i == n
            e += 1

            avg_l /= n

            println("y ~ $a * x + $b : ΔE[∇l] = $((avg_l-p_avg_l))")

            if e > 1 && (p_avg_l - avg_l) < 1e-3
                println("epochs = $e, E[∇l] = $avg_l")
                break
            end

            p_avg_l = avg_l
            avg_l = 0
        end

        i = mod1(i + 1, n)
    end

    lines(Ls)
end

