# Test time-domain state-space realization and simulation
@testset "State-space SISO" begin
    poles = [-5.0, -100 - 500im, -100 + 500im]
    residues = [2.0, 30 - 40im, 30 + 40im]
    d = 0.5
    D = d

    dt = 1e-5
    t_final = 1.0
    t = range(0, t_final, step = dt)  # time
    nt = length(t)
    input = ones(nt)
    input[1] = 0.0
    # Closed-formula step response
    step_response = ([
        sum(residues ./ poles .* (exp.(poles * t[i]) .- 1.0)) .+ d .* input[i] for
        i in eachindex(t)
    ])
    @assert isreal(step_response)
    step_response = real.(step_response)

    # Test Convolution
    yt = (convolution(poles, residues, input, dt, "recursive") .+ d .* input)
    @test isreal(yt)
    dif = (yt - step_response)
    @test norm(dif) < 0.022322
    @test maximum(abs.(dif)) < 0.0031768

    yt = (convolution(poles, residues, input, dt, "trapezoidal") .+ d .* input)
    @test isreal(yt)
    dif = (yt - step_response)
    @test norm(dif) < 0.0230
    @test maximum(abs.(dif)) < 0.0033

    # Real-Only
    A, B, C = rational_to_state_space(poles, residues; real_only = true, reduced = false)
    A_target = [
        -5.0 0.0 0.0
        0.0 -100.0 -500.0
        0.0 500.0 -100.0
    ]
    B_target = [1.0, 2.0, 0.0]
    C_target = [2.0 30.0 -40.0]
    @test A == A_target
    @test B == B_target
    @test C == C_target
    nx = size(A, 1)
    yt = simulate_state_space(A, B, C, D, 0, input, dt, nt)
    dif = (yt - step_response)
    @test norm(dif) < 0.00769
    @test maximum(abs.(dif)) < 0.00032

    # Complex
    A, B, C = rational_to_state_space(poles, residues; real_only = false, reduced = false)
    A_target = [
        -5.0+0.0im 0.0+0.0im 0.0+0.0im
        0.0+0.0im -100.0-500.0im 0.0+0.0im
        0.0+0.0im 0.0+0.0im -100.0+500.0im
    ]
    B_target = [1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im]
    C_target = [2.0 + 0.0im 30.0 - 40.0im 30.0 + 40.0im]
    @test A == A_target
    @test B == B_target
    @test C == C_target
    yt = simulate_state_space(A, B, C, D, 0, input, dt, nt)
    dif = (yt - step_response)
    @test norm(dif) < 0.00769
    @test maximum(abs.(dif)) < 0.00769

    # Reduced system
    A, B, C = rational_to_state_space(poles, residues; reduced = true)
    A_target = [
        -5.0+0.0im 0.0+0.0im
        0.0+0.0im -100.0+500.0im
    ]
    B_target = [1.0, 1.0]
    C_target = [2.0 + 0.0im 60.0 + 80.0im]
    @test A == A_target
    @test B == B_target
    @test C == C_target
    nx = size(A, 1)
    yt = simulate_state_space(A, B, C, D, 0, input, dt, nt)
    yt = real.(yt)  # discard imaginary part because a reduced model was used
    dif = (yt - step_response)
    @test norm(dif) < 0.00769
    @test maximum(abs.(dif)) < 0.00769
end
