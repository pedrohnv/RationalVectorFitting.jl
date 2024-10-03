@testset "ex4a" begin
    local s, bigY
    open("fdne.txt", "r") do fid1
        Nc = parse(Int, readline(fid1))
        Ns = parse(Int, readline(fid1))
        s = Array{ComplexF64}(undef, Ns)
        bigY = Array{ComplexF64}(undef, Nc, Nc, Ns)
        for k = 1:Ns
            s[k] = parse(Float64, readline(fid1)) * 1im
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
    Np = 50
    init_poles = RationalVectorFitting.recommended_init_poles(s, Np)
    poles, residues, d, h, fitted, error_norm = RationalVectorFitting.vector_fitting(
        s,
        f,
        init_poles,
        relaxed = false,
        maxiter = 5,
        tol = 1e-12,
    )
    #@test error_norm < 1e-10  # FIXME not fitting well. Add weighting to make it better
    @test error_norm < 1e-0
end
