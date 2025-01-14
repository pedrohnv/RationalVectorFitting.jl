#= Fitting of a 2x2 symmetric matrix representing a 2-inputs, 2-outputs system.

This is the example of the electrical circuit presented in section VIII of [1].

[1] B. Gustavsen and H. M. J. De Silva, "Inclusion of Rational Models in an
Electromagnetic Transients Program: Y-Parameters, Z-Parameters, S-Parameters,
Transfer Functions," in IEEE Transactions on Power Delivery, vol. 28, no. 2,
pp. 1164-1174, April 2013, doi: 10.1109/TPWRD.2013.2247067.
=#

using LinearAlgebra

@testset "matrix-fit" begin
    poles0 = [
        -0.476190502929121 + 0im,
        -79206.4441574531 + 0im,
        -1017.76384177315 - 3595.99255890622im,
        -1017.76384177315 + 3595.99255890622im,
        -875.957351183041 - 14593.7841813509im,
        -875.957351183041 + 14593.7841813509im,
        -1550.45082287519 - 34852.3036229938im,
        -1550.45082287519 + 34852.3036229938im,
        -15280.2844768572 - 122345.044171754im,
        -15280.2844768572 + 122345.044171754im,
    ]
    np = length(poles0)

    res11 = [
        0 + 0im,
        -198.785321394258 + 0im,
        0.219980455839476 + 0.215457612747253im,
        0.219980455839476 - 0.215457612747253im,
        81.6279911595347 + 6.47710195972654im,
        81.6279911595347 - 6.47710195972654im,
        333.260867183509 + 63.1092549627065im,
        333.260867183509 - 63.1092549627065im,
        4684.28382189821 + 239.040732309447im,
        4684.28382189821 - 239.040732309447im,
    ]
    res12 = [
        -2.15954408713931e-07 + 0im,
        -980.164975392196 + 0im,
        3.71765215165344 + 1.87852033203148im,
        3.71765215165344 - 1.87852033203148im,
        -105.402656273099 - 3.04541477612086im,
        -105.402656273099 + 3.04541477612086im,
        -111.883048304342 + 12.2553288771727im,
        -111.883048304342 - 12.2553288771727im,
        -3879.68279310344 - 1368.97367918464im,
        -3879.68279310344 + 1368.97367918464im,
    ]
    res22 = [
        47.6190547771491 + 0im,
        -4832.96941769953 + 0im,
        55.619118874847 + 9.01810901860201im,
        55.619118874847 - 9.01810901860201im,
        135.757522347641 - 2.90740586606814im,
        135.757522347641 + 2.90740586606814im,
        34.3218623196157 - 14.7282646981807im,
        34.3218623196157 + 14.7282646981807im,
        2921.31695569686 + 2118.58501483137im,
        2921.31695569686 - 2118.58501483137im,
    ]
    R = residues0 = zeros(ComplexF64, 2, 2, np)
    for i = 1:np
        R[1, 1, i] = res11[i]
        R[1, 2, i] = res12[i]
        R[2, 1, i] = res12[i]
        R[2, 2, i] = res22[i]
    end

    D0 = [
        0 0
        0 0.0833333333333335
    ]

    E0 = [
        0 0
        0 2e-7
    ]

    # Matrix Response F(s)
    ns = 500
    s = range(1, 2e5, length = ns) * 1im
    nc = 2
    Fs = zeros(ComplexF64, nc, nc, ns)
    for k = 1:ns
        for i2 = 1:nc
            for i1 = 1:nc
                R_sp = sum([residues0[i1, i2, p] / (s[k] - poles0[p]) for p = 1:np])
                Fs[i1, i2, k] = D0[i1, i2] + s[k] * E0[i1, i2] + R_sp
            end
        end
    end

    init_poles = recommended_init_poles(s, np)
    poles, residues, d, e, fitted, error_norm = symmetric_matrix_fitting(s, Fs, init_poles)

    @test error_norm < 8e-9
    @test norm(poles - poles0) < 1.6e-4
    @test norm(residues - residues0) < 8e-6
    @test norm(d - D0) < 3e-11
    @test norm(e - E0) < 2e-16
end
