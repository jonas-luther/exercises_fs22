
using LinearAlgebra

shift(s) = (s << 1 | (s >> 3 & 1)) & 0b1111

function T_mult(v)
    out = zeros(16)

    for j in 0:15
        out[shift(j) + 1] += v[j + 1]
    end
    
    out
end

function complex_to_string(z)
    if abs(z - 1.0) < 1e-6
        return ""
    end

    if abs(z + 1.0) < 1e-6
        return "-"
    end

    if abs(real(z)) < 1e-6 
        "$(imag(z))im"
    elseif abs(imag(z)) < 1e-6
        "$(real(z))"
    else
    "($z)"
    end
end

function state_to_string(s)
    str = ""

    nz = false

    for i in 1:16
        if abs(s[i]) > 1e-6

            coeff = complex_to_string(s[i])
            sstr = string([arr(c) for c in reverse(bitstring(i - 1))[1:4]]...)

            if nz
                if length(coeff) == 0 || coeff[1] != '-'
                    str *= " + $coeff|$sstr>"
                else
                    str *= " - $(coeff[2:end])|$sstr>"
                end
            else
                nz = true
                str *= "|$sstr>"
            end
        end
    end

    str
end

function orbit(s)
    orb = zeros(16, 4)

    orb[:, 1] .= s
    for i in 2:4
        orb[:, i] = T_mult(orb[:, i - 1])
    end

    orb[:, 1:rank(orb)]
end

function chi(orb)
    
    _, D = size(orb)

    chi = zeros(ComplexF64, 16, D)

    for i in 1:D

        k = (i - 1) * 4 / D

        p_k = 0.5 * (k) * pi

        for nu in 0:3
            chi[:, i] .+= exp(-im * p_k * nu) * orb[:, mod1(nu + 1, D)]
        end
    end

    chi
end

function basis_string(basis)
    str = ""
    r = rank(basis)
    str *= "{"
    for i in 1:r - 1
        str *= "$(state_to_string(basis[:,i])), "
    end
    str *= "$(state_to_string(basis[:,r]))}"
    str
end

function exercise_2()
    s = zeros(16)

    for i in 1:16
        s[i] = 1

        orb = orbit(s)
        println("orbit(|$(i - 1)>) = $(basis_string(orb))")
        chi_ = chi(orb)
        println("chi(|$(i - 1)>) = $(basis_string(chi_))")

        s[i] = 0
    end

   
end