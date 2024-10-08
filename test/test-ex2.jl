@testset "ex2" begin
    poles0 =
        2π .* [
            -41000,
            -4500,
            -100 - 5000im,
            -100 + 5000im,
            -120 - 15000im,
            -120 + 15000im,
            -3000 - 35000im,
            -3000 + 35000im,
            -200 - 45000im,
            -200 + 45000im,
            -1500 - 45000im,
            -1500 + 45000im,
            -500 - 70000im,
            -500 + 70000im,
            -1000 - 73000im,
            -1000 + 73000im,
            -2000 - 90000im,
            -2000 + 90000im,
        ]
    residues0 =
        2π .* [
            -83000,
            -3000,
            -5 - 7000im,
            -5 + 7000im,
            -20 - 18000im,
            -20 + 18000im,
            6000 - 45000im,
            6000 + 45000im,
            40 - 60000im,
            40 + 60000im,
            90 - 10000im,
            90 + 10000im,
            50000 - 80000im,
            50000 + 80000im,
            1000 - 45000im,
            1000 + 45000im,
            -5000 - 92000im,
            -5000 + 92000im,
        ]
    d0 = 0.2
    h0 = 2e-5
    freq = (range(0, 1e5, length = 200))
    s = 2im * pi * freq
    f = RationalVectorFitting.rational(s, poles0, residues0, d0, h0)
    init_poles =
        2π .* [
            -1e-2 + 1im,
            -1.11e2 + 1.11e4im,
            -2.22e2 + 2.22e4im,
            -3.33e2 + 3.33e4im,
            -4.44e2 + 4.44e4im,
            -5.55e2 + 5.55e4im,
            -6.66e2 + 7.77e4im,
            -8.88e2 + 8.88e4im,
            -1e3 + 1e5im,
        ]
    real_poles = filter(isreal, init_poles)
    complex_poles = filter(!isreal, init_poles)
    init_poles = sort!(
        [real_poles; complex_poles; conj(complex_poles)],
        by = RationalVectorFitting.cplxpair,
    )

    # No weighting
    poles, residues, d, h, fitted, error_norm =
        RationalVectorFitting.vector_fitting(s, f, init_poles; maxiter = 3, relaxed = false)
    @test error_norm < 1e-10

    poles, residues, d, h, fitted, error_norm =
        RationalVectorFitting.vector_fitting(s, f, init_poles; maxiter = 3, relaxed = true)
    @test error_norm < 1e-10

    # With weighting
    weight = @. 1.0 / sqrt(abs(f))
    poles, residues, d, h, fitted, error_norm = RationalVectorFitting.vector_fitting(
        s,
        f,
        init_poles,
        weight;
        maxiter = 3,
        relaxed = false,
    )
    @test error_norm < 1e-10

    poles, residues, d, h, fitted, error_norm = RationalVectorFitting.vector_fitting(
        s,
        f,
        init_poles,
        weight;
        maxiter = 3,
        relaxed = true,
    )
    @test error_norm < 1e-10
end
