using VectorFitting
using Test

begin  # ex3
    fid1 = split.(readlines("03PK10.txt"))
    f = zeros(ComplexF64, 160)
    for k = 1:160
        A1 = parse(Float64, fid1[k][1])
        A2 = parse(Float64, fid1[k][2])
        f[k] = A1 * exp(1im * A2 * pi / 180)
    end
    w = 2 * pi * range(0, 10e6, length=401)
    w = w[2:161]
    s = 1im .* w
    
    N = 12  # Order of approximation 
    init_poles = recommended_init_poles(s, N)
    poles, residues, d, h, fitted, error_norm = vector_fitting(s, f, init_poles; relaxed=false)
    @test error_norm < 1e-10
    #=
    freq = imag(s) ./ 2pi .* 1e-6
    p1 = plot(freq, abs.(f), xlabel="Frequency [MHz]", label="original", color=:blue)
    plot!(freq, abs.(fitted), label="fitted", color=:orange)
    plot!(freq, abs.(fitted .- f), label="deviation", yaxis=:log, legend=:left, color=:green)
    display(p1)

    p2 = plot(freq, angle.(f), xlabel="Frequency [MHz]", label="original")
    plot!(freq, angle.(fitted), label="fitted")
    display(p2)
    =#
end
