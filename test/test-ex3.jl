@testset "ex3" begin
    fid1 = split.(readlines("test/03PK10.txt"))
    f = zeros(ComplexF64, 160)
    for k = 1:160
        A1 = parse(Float64, fid1[k+1][1])
        A2 = parse(Float64, fid1[k+1][2])
        f[k] = A1 * exp(1im * A2 * pi / 180)
    end
    w = 2 * pi * range(0, 10e6, length = 401)
    w = w[2:161]
    s = 1im .* w
    N = 6  # Order of approximation
    init_poles = RationalVectorFitting.recommended_init_poles(s, N)

    # No weighting
    poles, residues, d, h, fitted, error_norm =
        RationalVectorFitting.vector_fitting(s, f, init_poles; maxiter = 5, relaxed = false)
    @test error_norm < 1e-1

    poles, residues, d, h, fitted, error_norm =
        RationalVectorFitting.vector_fitting(s, f, init_poles; maxiter = 5, relaxed = true)
    @test error_norm < 1e-1

    # With weighting
    weight = @. 1.0 / sqrt(abs(f))
    poles, residues, d, h, fitted, error_norm = RationalVectorFitting.vector_fitting(
        s,
        f,
        init_poles,
        weight;
        maxiter = 5,
        relaxed = false,
    )
    @test error_norm < 1e-1

    poles, residues, d, h, fitted, error_norm = RationalVectorFitting.vector_fitting(
        s,
        f,
        init_poles,
        weight;
        maxiter = 5,
        relaxed = true,
    )
    @test error_norm < 1e-1
end
