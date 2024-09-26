#= Fast Relaxed Vector Fitting

References
----------

[1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain
responses by vector fitting," in IEEE Transactions on Power Delivery, vol. 14,
no. 3, pp. 1052-1061, July 1999, doi: 10.1109/61.772353.

[2] B. Gustavsen, "Improving the pole relocating properties of vector fitting,"
in IEEE Transactions on Power Delivery, vol. 21, no. 3, pp. 1587-1592,
July 2006, doi: 10.1109/TPWRD.2005.860281.

[3] D. Deschrijver, M. Mrozowski, T. Dhaene and D. De Zutter, "Macromodeling of
Multiport Systems Using a Fast Implementation of the Vector Fitting Method,"
in IEEE Microwave and Wireless Components Letters, vol. 18, no. 6, pp. 383-385,
June 2008, doi: 10.1109/LMWC.2008.922585
=#
using LinearAlgebra
using Plots


"""To be used to sort an array by real, then complex conjugate pairs."""
function cplxpair(x)
    return (isreal(x), abs(imag(x)), real(x), imag(x))
end


function rational(s, poles, residues, d, h)
    return [sum(residues ./ (sk .- poles)) + d + sk * h for sk in s]
end


"""Builds a vector of recommended initial poles sorted by cplxpair."""
function recommended_init_poles(s, Npairs)
    s0 = imag(s[1])
    if isapprox(s0, 0.0)
        s0 = imag(s[2])
    end
    s1 = imag(s[end])
    init_poles = [(-0.01 + 1.0im) * sk for sk in range(s0, s1, length=Npairs÷2)]
    init_poles = sort!([init_poles; conj.(init_poles)], by=cplxpair)
    return init_poles
end


"""Builds the submatrix with the `1 / (s - p)`, `1.0` and `s` coefficients.
It is assumed that the poles are sorted by cplxpair."""
function build_subA!(A1, s, poles)
    Ns = length(s)
    Np = length(poles)
    skip_next = false
    for (i, p) in enumerate(poles)
        if skip_next
            skip_next = false
            continue
        elseif isreal(p)
            skip_next = false
            A1[1:Ns, i] .= 1.0 ./ (s .- p)
        else
            skip_next = true
            A1[1:Ns, i] .= 1.0 ./ (s .- p) + 1.0 ./ (s .- conj(p))
            A1[1:Ns, i+1] .= 1.0im ./ (s .- p) - 1.0im ./ (s .- conj(p))
        end
    end
    A1[1:Ns, Np+1] .= 1.0
    A1[1:Ns, Np+2] .= s
end


"""Stage 1 of the Vector Fitting."""
function pole_identification(s, f, poles, relaxed)
    Ns = length(s)
    Np = length(poles)
    Nres = Np + relaxed
    Ncols = Np + 2 + Nres
    A1_cplx = Array{ComplexF64}(undef, Ns, Ncols)
    Nrows = 2 * Ns + relaxed
    A1_reim = Array{Float64}(undef, Nrows, Ncols)
    Nc = (ndims(f) == 1) ? 1 : size(f)[2]
    A_sys = Array{Float64}(undef, (Nc * Nres), Nres)
    b_sys = zeros(Nc * Nres)
    @inline build_subA!(A1_cplx, s, poles)  # left block
    for n = 1:Nc
        A1_cplx[1:Ns, (Np + 3):Ncols] .= -f[1:Ns, n] .* A1_cplx[1:Ns, 1:Nres]  # right block
        A1_reim[1:Ns, :] .= real(A1_cplx)
        A1_reim[(Ns + 1):(2Ns), :] .= imag(A1_cplx)
        if relaxed && n == Nc
            A1_reim[end, 1:(Np + 2)] .= 0.0
            for i = 1:Nres
                A1_reim[end, Np + 2 + i] = real(sum(A1_cplx[:, i]))
            end
        end
        # Fast VF is a block-wise QR as we only want the last Nres values of
        # the solution. See [3].
        Q, R = qr!(A1_reim)
        i1 = (Np + 3)
        i2 = i1 + Np - 1 + relaxed
        k1 = 1 + (n - 1) * Nres
        k2 = k1 + Np - 1 + relaxed
        A_sys[k1:k2, :] .= R[i1:i2, i1:i2]
        if relaxed && n == Nc
            b_sys[k1:k2] .= Q[end, i1:i2] * Ns
        elseif !relaxed
            b_sys[k1:k2] .= transpose(Q[:, i1:i2]) * [real(f[1:Ns, n]); imag(f[1:Ns, n])]
        end
    end
    ldiv!(qr!(A_sys), b_sys)  # b = A \ b

    if relaxed
        sig_d = abs(b_sys[end])
        if sig_d < 1e-12
            b_sys[end] = 1e-8 * b_sys[end] / sig_d
            @warn "`d` of sigma too small at `iter. = $(iter)`. Consider stopping execution and setting `relaxed=false`. Resuming..."
        end
        b_sys[1:(end-1)] ./= b_sys[end]  # scale sigma's residues by its `d`
    end

    H = zeros(Np, Np)
    skip_next = false
    for (i, p) in enumerate(poles)
        if skip_next
            skip_next = false
            continue
        elseif isreal(p)
            skip_next = false
            H[:, i] .= -b_sys[i]
            H[i, i] += p
        else
            skip_next = true
            H[1:2:end, i] .= -2.0 * b_sys[i]
            H[1:2:end, i+1] .= -2.0 * b_sys[i+1]
            H[i, i] += real(p)
            H[i+1, i] += -imag(p)
            H[i, i+1] += imag(p)
            H[i+1, i+1] += real(p)
        end
    end
    return eigvals(H)
end


"""Stage 2 of the Vector Fitting."""
function residue_identification(s, f, poles)
    Ns = length(s)
    Np = length(poles)
    Nc = (ndims(f) == 1) ? 1 : size(f)[2]
    residues = Array{ComplexF64}(undef, Np, Nc)
    d = zeros(Nc)
    h = similar(d)
    Nrows = 2 * Ns
    Ncols = Np + 2
    A1_cplx = Array{ComplexF64}(undef, Ns, Ncols)
    A_sys = Array{Float64}(undef, Nrows, Ncols)
    X_sys = Array{Float64}(undef, Ncols, Nc)

    @inline build_subA!(A1_cplx, s, poles)
    A_sys[1:Ns, :] .= real(A1_cplx)
    A_sys[(Ns+1):end, :] .= imag(A1_cplx)
    X_sys = A_sys \ [real(f); imag(f)]
    for n = 1:Nc
        skip_next = false
        for (i, p) in enumerate(poles)
            if skip_next
                skip_next = false
                continue
            elseif isreal(p)
                skip_next = false
                residues[i, n] = X_sys[i, n]
            else
                skip_next = true
                residues[i, n] = complex(X_sys[i, n], X_sys[i + 1, n])
                residues[i + 1, n] = conj(residues[i, n])
            end
        end
        d[n] = X_sys[Np+1, n]
        h[n] = X_sys[Np+2, n]
    end
    return residues, d, h
end



"""Perform the vector fitting of the array `f` with complex frequency `s`
using a set of initial poles. It is assumed that `init_poles`.

`f` can be a matrix of dimensions `(Ns, Nc)` and the fitting will be over
its columns with a set of common poles.

`relaxed` controls the nontriviality constraint. See [2].

`force_stable` controls if unstable poles should be reflected to the semi-left
complex plane.
"""
function vector_fitting(s, f, init_poles; relaxed=true, force_stable=true, maxiter=20, tol=1e-12)
    if !allequal(real(s))
        throw(error("It is expected that `allequal(real(s)) == true`"))
    end

    if any(imag(s) .< 0.0)
        throw(error("It is expected that `all(imag(s) .>= 0) == true`"))
    end
    
    if ndims(f) == 1
        Nc = 1
    elseif ndims(f) == 2
        Nc = size(f)[2]
    else
        throw(error("It is expected `f` to have 1 or 2 dimensions."))
    end

    Ns = length(s)
    if Ns != size(f)[1]
        throw(error("`f` must have the same number of rows as `s`"))
    end

    poles = sort!(complex(init_poles), by=cplxpair)
    fitted = similar(f)
    error_norm = Inf
    local residues, d, h
    for iter = 1:maxiter
        if error_norm < tol
            println("convergence achieved at iter. = $(iter)")
            println("error_norm = $(error_norm)")
            break
        end
        poles = pole_identification(s, f, poles, relaxed)
        if force_stable
            for (i, p) in enumerate(poles)
                re_p, im_p = reim(p)
                if re_p > 0.0
                    poles[i] = complex(-re_p, im_p)
                end
            end
        end
        residues, d, h = residue_identification(s, f, poles)
        for n = 1:Nc
            fitted[:,n] .= rational(s, poles, residues[:,n], d[n], h[n])
        end
        error_norm = norm(f .- fitted, 2)
    end
    perm = sortperm(poles, by=cplxpair)
    poles = poles[perm]
    for n = 1:Nc
        residues[:, n] = residues[perm, n]
    end
    return poles, residues, d, h, fitted
end


begin  # ex1
    Ns = 101
    s = 2im * pi * exp10.(range(0, 4, length=Ns))
    poles0 = [-5.0, -100 - 500im, -100 + 500im]
    residues0 = [2.0, 30 - 40im, 30 + 40im]
    d0 = 0.5
    h0 = 0.0
    f = rational(s, poles0, residues0, d0, h0)
    init_poles = -2pi * exp10.(range(0, 4, length=3))
    poles, residues, d, h, fitted = vector_fitting(s, f, init_poles)
    freq = imag(s) ./ (2pi) * 1e-3
    p1 = plot(freq, abs.(f), xlabel="Frequency [kHz]", label="original", color=:blue)
    plot!(freq, abs.(fitted), label="fitted", color=:orange, linestyle=:dash)
    plot!(freq, abs.(fitted .- f), label="deviation", xaxis=:log, yaxis=:log, legend=:left, color=:green)
    display(p1)
end



begin  # ex2
    poles0 = 2π .* [
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
    residues0 = 2π .* [
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

    freq = (range(0, 1e5, length=200))
    s = 2im * pi * freq
    f = rational(s, poles0, residues0, d0, h0)
    # insert noise:
    #using Random, Distributions
    #f .+= rand(Normal(0.0, sqrt(10)), length(s))

    init_poles = 2π .* [
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
    init_poles = sort!([real_poles; complex_poles; conj(complex_poles)], by=cplxpair)
    #init_poles = recommended_init_poles(s, 20)
    poles, residues, d, h, fitted = vector_fitting(s, f, init_poles)

    p1 = plot(freq .* 1e-3, abs.(f), xlabel="Frequency [kHz]", label="original")
    plot!(freq .* 1e-3, abs.(fitted), label="fitted")
    plot!(freq .* 1e-3, abs.(fitted .- f), label="deviation", yaxis=:log, legend=:left)
    display(p1)
end



begin  # ex3
    fid1 = split.(readlines("03PK10.txt"))
    f = zeros(ComplexF64, 160)
    for k = 1:160
        A1 = parse(Float64, fid1[k][1])
        A2 = parse(Float64, fid1[k][2])
        f[k] = A1 * exp(1im * A2 * pi / 180)
    end
    w = 2 * pi * range(0, 10e6, length=401)
    w = w[2:161]
    s = 1im .* w
    
    N = 12  # Order of approximation 
    init_poles = recommended_init_poles(s, N)
    poles, residues, d, h, fitted = vector_fitting(s, f, init_poles; relaxed=false)
    
    freq = imag(s) ./ 2pi .* 1e-6
    p1 = plot(freq, abs.(f), xlabel="Frequency [MHz]", label="original", color=:blue)
    plot!(freq, abs.(fitted), label="fitted", color=:orange)
    plot!(freq, abs.(fitted .- f), label="deviation", yaxis=:log, legend=:left, color=:green)
    display(p1)

    p2 = plot(freq, angle.(f), xlabel="Frequency [MHz]", label="original")
    plot!(freq, angle.(fitted), label="fitted")
    display(p2)
end



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
                    A1 = parse(Float64, readline(fid1))
                    a2 = parse(Float64, readline(fid1))
                    bigY[row,col,k] = A1 + a2 * 1im
                end
            end
        end
    end
    f = transpose(bigY[:,1,:])
    Ns, Nc = size(f)
    Np = 50
    init_poles = recommended_init_poles(s, Np)
    poles, residues, d, h, fitted = vector_fitting(s, f, init_poles, relaxed=false, maxiter=5, tol=1e-12)
    @show norm(f .- fitted, 2)

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
end

# TODO inserir pesos
