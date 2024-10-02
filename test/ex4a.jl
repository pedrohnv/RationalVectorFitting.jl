using VectorFitting
using Test

begin  # ex4a
    open("fdne.txt", "r") do fid1
        Nc = parse(Int, readline(fid1))
        Ns = parse(Int, readline(fid1))
        global s = Array{ComplexF64}(undef, Ns)
        global bigY = Array{ComplexF64}(undef, Nc, Nc, Ns)
        for k = 1:Ns
            s[k] = parse(Float64, readline(fid1)) * 1im
            for row = 1:Nc
                for col = 1:Nc
                    a1 = parse(Float64, readline(fid1))
                    a2 = parse(Float64, readline(fid1))
                    bigY[row,col,k] = complex(a1, a2)
                end
            end
        end
    end
    f = transpose(bigY[:,1,:])
    Ns, Nc = size(f)
    Np = 50
    init_poles = recommended_init_poles(s, Np)
    poles, residues, d, h, fitted, error_norm = vector_fitting(s, f, init_poles, relaxed=false, maxiter=5, tol=1e-12)
    @test error_norm < 1e-10
    #=
    freq = imag(s) ./ 2pi

    p1 = plot(freq, abs.(f), xlabel="Frequency [Hz]", label="", color=:blue)
    plot!(freq, abs.(fitted), label="", color=:orange)#, linestyle=:dash)
    plot!(freq, abs.(fitted .- f), label="", color=:green, yaxis=:log, legend=:topleft)
    plot!(freq[1:2], abs.(f)[1:2], label="original", color=:blue)
    plot!(freq[1:2], abs.(fitted)[1:2], label="fitted", color=:orange)
    plot!(freq[1:2], abs.(fitted .- f)[1:2], label="deviation", color=:green)
    display(p1)

    using DSP
    function continuos_angle(A)
        return mapslices(unwrap, angle.(A), dims=1)
    end
    p2 = plot(freq, continuos_angle(f), xlabel="Frequency [Hz]", ylabel="Phase angle [rad]", label="", color=:blue)
    plot!(freq, continuos_angle(fitted), label="", color=:orange)
    plot!(freq[1:2], angle.(f)[1:2], label="original", color=:blue)
    plot!(freq[1:2], angle.(fitted)[1:2], label="fitted", color=:orange)
    display(p2)
    =#
end
