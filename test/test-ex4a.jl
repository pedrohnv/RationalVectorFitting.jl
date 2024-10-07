@testset "ex4a" begin
    local s, bigY
    open("test/fdne.txt", "r") do fid1
        Nc = parse(Int, readline(fid1))
        Ns = parse(Int, readline(fid1))
        s = Array{ComplexF64}(undef, Ns)
        bigY = Array{ComplexF64}(undef, Nc, Nc, Ns)
        for k = 1:Ns
            s[k] = complex(0, parse(Float64, readline(fid1)))
            for row = 1:Nc
                for col = 1:Nc
                    a1 = parse(Float64, readline(fid1))
                    a2 = parse(Float64, readline(fid1))
                    bigY[row, col, k] = complex(a1, a2)
                end
            end
        end
    end
    f = transpose(bigY[:, 1, :])
    Ns, Nc = size(f)
    Np = 100  # FIXME 50 should be enough, but is giving a bad fit...
    init_poles = RationalVectorFitting.recommended_init_poles(s, Np)
    weight = @. 1.0 / sqrt(abs(f))

    poles, residues, d, h, fitted, error_norm = RationalVectorFitting.vector_fitting(
        s,
        f,
        init_poles,
        weight;
        maxiter = 6,
        relaxed = false,
        tol = 0.1,
    )

    init_poles = poles
    poles, residues, d, h, fitted, error_norm = RationalVectorFitting.vector_fitting(
        s,
        f,
        init_poles,
        weight;
        maxiter = 6,
        relaxed = true,
        tol = 0.1,
    )

    @test error_norm < 1e-1
end
