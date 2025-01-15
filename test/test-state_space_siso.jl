# Test time-domain state-space realization and simulation
@testset "State-space SISO" begin
    poles = [-5.0, -100 - 500im, -100 + 500im]
    residues = [2.0, 30 - 40im, 30 + 40im]
    D = 0.5

    dt = 1e-4
    t_final = 1.0
    t = range(0, t_final, step = dt)  # time
    nt = length(t)
    input = ones(nt)
    input[1] = 0.0
    # Closed-formula step response
    step_response = ([
        sum(residues ./ poles .* (exp.(poles * t[i]) .- 1.0)) .+ D .* input[i] for
        i in eachindex(t)
    ])
    @assert isreal(step_response)
    step_response = real.(step_response)

    # Test Convolution
    yt = (convolution(poles, residues, input, dt, "recursive") .+ D .* input)
    @test isreal(yt)
    dif = (yt - step_response)
    @test norm(dif) < 0.022322
    @test maximum(abs.(dif)) < 0.0031768

    yt = (convolution(poles, residues, input, dt, "trapezoidal") .+ D .* input)
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
    x0 = zeros(nx)
    yt = simulate_state_space(A, B, C, D, input, x0, dt, nt)
    dif = (yt - step_response)
    @test norm(dif) < 0.000769
    @test maximum(abs.(dif)) < 7.97e-5

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
    yt = simulate_state_space(A, B, C, D, input, x0, dt, nt)
    dif = (yt - step_response)
    @test norm(dif) < 0.000769
    @test maximum(abs.(dif)) < 0.000769

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
    x0 = zeros(nx)
    yt = simulate_state_space(A, B, C, D, input, x0, dt, nt)
    yt = real.(yt)  # discard imaginary part because a reduced model was used
    dif = (yt - step_response)
    @test norm(dif) < 0.000769
    @test maximum(abs.(dif)) < 0.000769
end
