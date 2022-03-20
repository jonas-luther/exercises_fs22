module QM

using LinearAlgebra


export s_x, s_y, s_z

const s_x = [0 1; 
             1 0]
const s_y = [0 1im; 
             -1im 0]
const s_z = [1 0;
             0 -1]

function plot_wave_function(xs, psi)

    lines(xs, real.(psi), label="Re ψ")
    lines!(xs, imag.(psi), label="Im ψ")
    ylims!(-2, 2)
    axislegend()
    current_figure()

end

end