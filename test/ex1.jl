using VectorFitting
using Test

begin  # ex1
    Ns = 101
    s = 2im * pi * exp10.(range(0, 4, length=Ns))
    poles0 = [-5.0, -100 - 500im, -100 + 500im]
    residues0 = [2.0, 30 - 40im, 30 + 40im]
    d0 = 0.5
    h0 = 0.0
    f = rational(s, poles0, residues0, d0, h0)
    init_poles = -2pi * exp10.(range(0, 4, length=3))
    poles, residues, d, h, fitted, error_norm = vector_fitting(s, f, init_poles)
    @test error_norm < 1e-10
    #=
    freq = imag(s) ./ (2pi) * 1e-3
    p1 = plot(freq, abs.(f), xlabel="Frequency [kHz]", label="original", color=:blue)
    plot!(freq, abs.(fitted), label="fitted", color=:orange, linestyle=:dash)
    plot!(freq, abs.(fitted .- f), label="deviation", xaxis=:log, yaxis=:log, legend=:left, color=:green)
    display(p1)
    =#
end