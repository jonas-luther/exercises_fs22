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

    vals = Vector{Float64}()

    max_size = 0
    total = 0

    for k in 0:N - 1
        lam = exp(2π * im * k / N)

        orbs = Vector{Int64}()
        r2o = Dict{Int64,Int64}()

        for (repr, orb) in orbits
            if abs(lam^length(orb) - 1) < 1e-6
                push!(orbs, repr)
                r2o[repr] = length(orbs)
            end
        end
        
        size = length(orbs)
        H = zeros(ComplexF64, size, size)

        for i in 1:size
            s = orbs[i]
            H[i,i] += (2 * count_ones(~(s ⊻ (s >> 1 | (s & 1) << (N - 1))) & (2^N - 1)) - N)
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
                        beta = r2o[gamma]
                        n_beta = length(orbits[s_beta])
 
                        H[beta, alpha] += lam^(j - zeta) / sqrt(n_alpha * n_beta)
                    end
                end
            end
        end
        
        max_size = max(size, max_size)
        total += size

        append!(vals, eigvals(Hermitian(H)))
    end

    sort!(vals)

    dedup = [vals[1]]

    for i in 2:n_states
        if abs(vals[i] - dedup[end]) > 1e-6
            push!(dedup, vals[i])
        end
    end

    dedup

end

function test_block(N)
    @time m = main(N)

    scatter(m, label="Block Diagonal")

    axislegend()

    current_figure()
end


