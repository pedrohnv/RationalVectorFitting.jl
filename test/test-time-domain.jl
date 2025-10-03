# Test the Time-Domain Vector Fitting
using Test
using LinearAlgebra
using RationalVectorFitting.TimeDomain
#=
@testset "TimeDomainVF Basic Functionality" begin

    @testset "rational_to_state_space - Real Poles" begin
        poles = [-1.0, -2.0, -3.0]
        residues = [1.0, 2.0, 3.0]

        A, B, C = TimeDomain.rational_to_state_space(poles, residues)

        @test size(A) == (3, 3)
        @test size(B) == (3,)
        @test size(C) == (1, 3)
        @test all(isreal, A) && all(isreal, B) && all(isreal, C)
        @test diag(A) ≈ poles
        @test B ≈ [1.0, 1.0, 1.0]
        @test C ≈ [1.0 2.0 3.0]
    end

    @testset "rational_to_state_space - Complex Poles" begin
        poles = [-1+2im, -1-2im, -3.0]
        residues = [1+0.5im, 1-0.5im, 2.0]

        A, B, C = TimeDomain.rational_to_state_space(poles, residues; real_only = true)

        @test size(A) == (3, 3)
        @test all(isreal, A) && all(isreal, B) && all(isreal, C)
        # Test that complex poles are properly converted to real 2x2 blocks
        @test A[1, 1] ≈ -1.0
        @test A[1, 2] ≈ 2.0
        @test A[2, 1] ≈ -2.0
        @test A[2, 2] ≈ -1.0
        @test A[3, 3] ≈ -3.0
    end

    @testset "rational_to_state_space - Matrix Residues" begin
        poles = [-1.0, -2.0]
        residues = [1.0 2.0; 3.0 4.0]  # 2 outputs, 2 poles

        A, B, C = TimeDomain.rational_to_state_space(poles, residues)

        @test size(A) == (2, 2)
        @test size(B) == (2,)
        @test size(C) == (2, 2)
        @test diag(A) ≈ poles
        @test C ≈ residues
    end

    @testset "rational_to_state_space - Reduced Option" begin
        poles = [-1+2im, -1-2im, -3.0]
        residues = [1+0.5im, 1-0.5im, 2.0]

        A, B, C = TimeDomain.rational_to_state_space(poles, residues; reduced = true)

        # Should only keep poles with non-negative imaginary part
        @test size(A) == (2, 2)  # One complex pole + one real pole
        @test A[1, 1] ≈ -1+2im
        @test A[2, 2] ≈ -3.0
        # Residue for complex pole should be doubled
        @test C[1] ≈ 2.0 + 1.0im  # 2 * (1+0.5im)
    end
end

@testset "symmetric_rational_to_state_space" begin
    @testset "Symmetric 2x2 System" begin
        poles = [-1.0, -2.0]
        # Symmetric residues for 2x2 system
        residues = cat([1.0 0.5; 0.5 2.0], [2.0 1.0; 1.0 3.0], dims = 3)

        A, B, C = TimeDomain.symmetric_rational_to_state_space(poles, residues)

        @test size(A) == (4, 4)  # 2 poles × 2 states
        @test size(B) == (4, 2)  # 2 inputs
        @test size(C) == (2, 4)  # 2 outputs
    end

    @testset "Symmetric Validation" begin
        poles = [-1.0]
        residues = cat([1.0 2.0; 3.0 4.0], dims = 3)  # Not symmetric!

        @test_throws ErrorException TimeDomain.symmetric_rational_to_state_space(
            poles,
            residues,
        )
    end
end

@testset "simulate_state_space" begin
    @testset "Simple RC Circuit Analog" begin
        # Simple first-order system: dy/dt = -y + u
        # with u(t) = 1, then y(t) = 1 - exp(-t)
        A = [-1.0]
        B = [1.0]
        C = [1.0]
        D = [0.0]
        E = [0.0]

        dt = 0.001
        nt = 1000
        t = range(0, dt*(nt-1), length = nt)
        u = ones(nt)  # Step input
        y = TimeDomain.simulate_state_space(A, B, C, D, E, u, dt, nt)

        @test size(y) == (nt, 1)
        analytical = t -> 1 - exp(-t)
        err = @. abs(y - analytical(t + dt/2))
        @test maximum(err) < 1e-6
    end

    @testset "Impulse Response" begin
        # System with direct feedthrough
        # dy/dt = -2y + u
        # with u(t) = δ(t), then
        # y(t) = 0.5 * δ(t) + exp(-2t)
        A = [-2.0]
        B = [1.0]
        C = [1.0]
        D = [0.5]
        E = [0.0]

        dt = 0.001
        nt = 100
        t = range(0, dt*(nt-1), length = nt)
        u = zeros(nt)
        x0 = B  # For impulse response, initial condition x(0⁺) = B

        y = TimeDomain.simulate_state_space(A, B, C, D, E, u, dt, nt, x0)

        @test y[1] ≈ 1.0 atol=1e-9  # Immediate response from D
        analytical = t -> exp(-2t)
        err = @. abs(y - analytical(t + dt/2))
        @test maximum(err) < 1e-3
    end
end

@testset "convolution" begin
    @testset "Exponential Decay" begin
        poles = [-1.0]
        residues = [1.0]
        dt = 0.1
        nt = 20
        t = range(0, dt*(nt-1), length = nt)
        yt = ones(nt)  # Step function

        result = TimeDomain.convolution(dt, yt, poles, residues, "recursive")

        @test length(result) == nt
        # Convolution of step with exp(-t) should give (1 - exp(-t))
        expected = 1.0 .- exp.(-t)
        @test result ≈ expected atol=0.1
    end

    @testset "Multiple Poles" begin
        poles = [-1.0, -2.0]
        residues = [1.0, 2.0]
        dt = 0.1
        nt = 10
        yt = ones(nt)

        result = TimeDomain.convolution(dt, yt, poles, residues, "trapezoidal")

        @test length(result) == nt
        # Should be smooth and increasing
        @test all(diff(real(result)) .> 0)
    end

    @testset "Invalid Formula" begin
        poles = [-1.0]
        residues = [1.0]
        yt = [1.0]
        dt = 0.1

        @test_throws ErrorException TimeDomain.convolution(
            dt,
            yt,
            poles,
            residues,
            "invalid_formula",
        )
    end
end

@testset "Error Handling" begin
    @testset "Dimension Mismatch" begin
        poles = [-1.0, -2.0]
        residues = [1.0]  # Wrong size

        @test_throws ErrorException TimeDomain.rational_to_state_space(poles, residues)
    end

    @testset "Invalid Residue Dimensions" begin
        poles = [-1.0]
        residues = cat([1 2 3], dims = 3)  # 3D but wrong

        @test_throws ErrorException TimeDomain.rational_to_state_space(poles, residues)
    end
end
=#


dt = 0.02
nt = 250
t = range(0, dt*(nt-1), length = nt)

# Known model
poles = [-2.0, -1.0 + 5.0im, -1.0 - 5.0im]
residues = [1.0, 0.5 + 0.2im, 0.5 - 0.2im]
d = 0.1
h = 0.0
A, B, C = rational_to_state_space(poles, residues; real_only = true)
D = [d]
E = [h]
u = ones(nt)  # Step input
y_known = simulate_state_space(A, B, C, D, E, u, dt, nt)
plot(t, y_known)


init_poles = [-10, -100, -1000.0]
vin = zeros(nt)
vin[2] = 1 / dt
vin .+= 1e-12
vout = y_known
res = TimeDomain.vector_fitting_time_domain(
    dt,
    vin,
    vout,
    init_poles,
    has_direct_feedthrough = true,
    niter = 1,
    formula = "recursive",
)
qpol, fitted, pointwise_rmsd, rmsd, pointwise_mean_abs_d, mean_abs_d = res
qpol


using Plots

fid1 = split.(readlines("test/td_surge.csv"))
N = length(fid1) - 1
t = zeros(N)
f = zeros(N)
for k = 1:N
    val = split(fid1[k+1][1], ",")
    t[k] = parse(Float64, val[1])
    f[k] = parse(Float64, val[2])
end
