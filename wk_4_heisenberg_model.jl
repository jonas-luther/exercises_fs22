using LinearAlgebra

function roll(N, n, b) 
    s = Int64(b)
    (s << n | (s >> (N - n))) & (1 << N - 1)
end  

function state_to_str(N, s)
    bitstring(s)[end - N + 1:end]
end

function main(N)
    n_states = 1 << N

    orep = Set{Int64}() # orbit representatives
    s2r = zeros(Int64, n_states) # maps state index to orbit representative index
    s2d = zeros(Int64, n_states) # maps state index to representative distance
    orbits = Dict{Int64,Set{Int64}}()

    for state in 0:(n_states - 1)
        new = true
        rep = -1
        for i in 0:N - 1
            rsn = roll(N, i, state)  

            if rsn in orep
                rep = rsn
                new = false
                s2d[state + 1] = N - i
                break
            end
        end

        if new
            push!(orep, state)
            rep = state
            orbits[rep] = Set()
            s2d[state + 1] = 0
        end

        push!(orbits[rep], state)
        s2r[state + 1] = rep
    end

    for (repr, s) in orbits
        print("($repr)")
        print(bitstring(repr)[end - N + 1:end])
        print(" : ")

        for state in s
            print(" $(state_to_str(N, state)) : ($(s2r[state + 1]), $(s2d[state + 1])),")
        end

        println()
    end


    for s in 0:n_states - 1
        s_r = roll(N, s2d[s + 1], s2r[s + 1])
        if s != s_r
            println("$(state_to_str(N, s)) != $(state_to_str(N, s_r))")
            println("rep = $(state_to_str(N, s2r[s + 1]))")
            println("distance = $(s2d[s + 1])")
        end
    end

    vals = Vector{ComplexF64}()

    max_size = 0
    total = 0

    for k in 0:N - 1
        lam = exp(2π * im * k / N)

        orbs = Vector{Int64}()

        println(k == 0 ? "λ = 1:" : "λ = exp(2πi * $k/$N):")

        for (repr, orb) in orbits
            if abs(lam^length(orb) - 1) < 1e-6
                push!(orbs, repr)
            end
        end
        
        size = length(orbs)
        H = zeros(ComplexF64, size, size)

        for i in 1:size
            s = orbs[i]
            H[i,i] += (2 * count_ones(~(s ⊻ (s >> 1 | (s & 1) << (N - 1))) & (2^N - 1)) - N)
            println("- $(bitstring(s)[end - N + 1:end])\n E_z = $(real(H[i,i]))")
        end
  
        for alpha in 1:size
            s_alpha = orbs[alpha]
            n_alpha = length(orbits[s_alpha])

            for v in 0:(N - 1)     
                for j in 0:n_alpha - 1
                    s_j = roll(N, j, s_alpha)

                    s_flip = s_j ⊻ (1 << v)

                    zeta = s2d[s_flip + 1]
                    gamma = s2r[s_flip + 1]

                    s_beta = gamma
                    
                    if s_beta in orbs
                        beta = indexin(s_beta, orbs)[1]
                        n_beta = length(orbits[s_beta])
 
                        H[beta, alpha] += lam^(j - zeta) / sqrt(n_alpha * n_beta)
                    end
                end
            end
        end
        
        max_size = max(size, max_size)
        total += size
        
        display(H)

        println(ishermitian(H))
        
        append!(vals, eigvals(H))
    end
    
    println("largest block size = $max_size vs $(2^N)")
    println("total = $total")
    
    
    sort(real.(vals))
    end

function compare(N)
    @time t = testDense(N)
    @time m = main(N)

    scatter(t, label="Dense")
    scatter!(m, label="Block Diagonal")

    axislegend()

    current_figure()
end


