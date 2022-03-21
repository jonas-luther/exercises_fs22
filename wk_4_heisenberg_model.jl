using LinearAlgebra
using GLMakie

function roll(N, n, b) 
    s = Int64(b)
    (s << n | (s >> (N-n))) & (1<<N - 1)
end  

function main(N)
    n_states = 1 << N

    orep = Set{Int64}() # orbit representatives
    s2r = zeros(Int64, n_states) # maps state index to orbit representative index
    s2d = zeros(Int64, n_states) # maps state index to representative distance
    orbits = Dict{Int64, Set{Int64}}()

    for state in 0:(n_states - 1)
        new = true
        rep = -1
        rsn = state
        for i in 0:N-1
            if rsn in orep
                rep = rsn
                new = false
                s2d[state] = i
                break
            end

            rsn = roll(N, 1, rsn)  
        end

        if new
            push!(orep, state)
            rep = state
            orbits[rep] = Set()
            s2d[state+1] = 0
        end
        push!(orbits[rep], state)
        s2r[state+1] = rep
    end

    for (repr, s) in orbits
        print(bitstring(repr)[end-N+1:end])
        print(" : ")

        for state in s
            print(bitstring(state)[end-N+1:end])
            print(", ")
        end

        println()
    end

    vals = Vector{Float64}()

    max_size = 0
    total = 0

    for k in 0:N-1
        lam = exp(2π * im * k / N)

        orbs = Vector{Int64}()

        println(k == 0 ? "λ = 1:" : "λ = exp(2πi * $k/$N):")

        for (repr, orb) in orbits
            if abs(lam^length(orb) - 1) < 1e-6
                push!(orbs, repr)
               
            end
        end
        
        size = length(orbs)
        H = zeros(size, size)

        for i in 1:size
            s = orbs[i]
            H[i,i] += (2 * count_ones(~(s ⊻ (s >> 1 | (s & 1) << (N - 1))) & (2^N - 1)) - N)
            println("- $(bitstring(s)[end-N+1:end])\n E_z = $(real(H[i,i]))")
        end

        # ??? confusion
        # TODO: do the thing
        # profit

        for beta in 1:size
            s_beta = orbs[beta]
            
            for v in 0:(N-1)
                for alpha in 1:size
                    
                    s_alpha = orbs[alpha]

                    s_flip = s_alpha ⊻ (1 << v)

                    zeta = s2d[s_flip+1]
                    gamma = s2r[s_flip+1]

                    for j in 0:size-1
                        if gamma == s_beta
                            H[beta, alpha] += lam^(j - zeta) / size
                        end
                    end
                end
            end

        end

        
        max_size = max(size, max_size)
        total += size

        display(H)

        println(ishermitian(H))
        
        
        union!(vals, eigvals(Hermitian(H)))
    end

    println("largest block size = $max_size vs $(2^N)")
    println("total = $total")

    println(vals)
    scatter(zeros(size(vals)), vals)
end




