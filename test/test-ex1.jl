# Uses one of the MATLAB examples as a test.
@testset "ex1" begin
    Ns = 101
    s = 2im * pi * exp10.(range(0, 4, length = Ns))
    poles0 = [-5.0, -100 - 500im, -100 + 500im]
    residues0 = [2.0, 30 - 40im, 30 + 40im]
    d0 = 0.5
    h0 = 0.0
    f = RationalVectorFitting.rational(s, poles0, residues0, d0, h0)
    init_poles = -2pi * exp10.(range(0, 4, length = 3))

    # No weighting
    poles, residues, d, h, fitted, error_norm =
        RationalVectorFitting.vector_fitting(s, f, init_poles; maxiter = 1, relaxed = false)
    @test error_norm < 1e-12

    poles, residues, d, h, fitted, error_norm =
        RationalVectorFitting.vector_fitting(s, f, init_poles; maxiter = 1, relaxed = true)
    @test error_norm < 1e-12

    # With weighting
    weight = @. 1.0 / sqrt(abs(f))
    poles, residues, d, h, fitted, error_norm = RationalVectorFitting.vector_fitting(
        s,
        f,
        init_poles,
        weight;
        maxiter = 1,
        relaxed = false,
    )
    @test error_norm < 1e-12

    poles, residues, d, h, fitted, error_norm = RationalVectorFitting.vector_fitting(
        s,
        f,
        init_poles,
        weight;
        maxiter = 1,
        relaxed = true,
    )
    @test error_norm < 1e-12
end
