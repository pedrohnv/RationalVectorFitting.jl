module RationalVectorFitting

export rational,
    cplxpair,
    recommended_init_poles,
    pole_identification,
    residue_identification,
    vector_fitting


using LinearAlgebra


"""
    cplxpair(x)

To be used to sort an array by real values, then complex conjugate pairs.
The more positive reals will be first, then the pairs with smaller imaginary part.

# Example
```julia
real_values = [-1, -3, -2]
complex_values = [-1 + 1im, -2 - 2im, -1 + 2im]
sorted_values = sort!([complex_values; conj(complex_values); real_values], by = cplxpair)

# output

9-element Vector{Complex{Int64}}:
 -1 + 0im
 -2 + 0im
 -3 + 0im
 -1 - 1im
 -1 + 1im
 -1 - 2im
 -1 + 2im
 -2 - 2im
 -2 + 2im
```
"""
function cplxpair(x)
    return (!isreal(x), abs(imag(x)), abs(real(x)), imag(x))
end


@doc raw"""
    rational(s, poles, residues, d, h) -> `f(s)`

Rational transfer function with complex frequencies `s`, a set of poles `a_n`,
residues `r_n` and real constants `d` and `h`.

```math
\sum_{n=1}^N \frac{r_n}{s - a_n} + d + s h
```
"""
function rational(s, poles, residues, d, h)
    return [sum(residues ./ (sk .- poles)) + d + sk * h for sk in s]
end


"""
    recommended_init_poles(s, Npairs) -> init_poles

Builds a vector of recommended initial poles sorted by [`cplxpair`](@ref).
"""
function recommended_init_poles(s, Npairs)
    s0 = imag(s[1])
    if isapprox(s0, 0.0)
        s0 = imag(s[2])
    end
    s1 = imag(s[end])
    init_poles = [(-0.01 + 1.0im) * sk for sk in range(s0, s1, length = Npairs ÷ 2)]
    init_poles = sort!([init_poles; conj.(init_poles)], by = cplxpair)
    return init_poles
end


"""
    build_Abase!(A1, s, poles)

Builds the base matrix with the `1 / (s - p)`, `1.0` and `s` coefficients.
It is assumed that the poles were sorted by [`cplxpair`](@ref).
"""
@inline function build_Abase!(A1, s, poles)
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


"""
    pole_identification(s, f, poles, weight, relaxed) -> new_poles

Stage 1 of the Vector Fitting.

See also [`vector_fitting`](@ref), [`residue_identification`](@ref).
"""
function pole_identification(s, f, poles, weight, relaxed)
    Ns = length(s)
    Np = length(poles)
    Nres = Np + relaxed
    Ncols = Np + 2 + Nres
    Nrows = 2 * Ns + relaxed
    A1_reim = zeros(Nrows, Ncols)
    Nc = (ndims(f) == 1) ? 1 : size(f)[2]
    A_sys = zeros((Nc * Nres), Nres)
    b_sys = zeros(Nc * Nres)
    x_sys = view(b_sys, 1:Nres)
    scale = sqrt(sum([norm(weight[1:Ns, n] .* f[1:Ns, n])^2 for n = 1:Nc]))
    scale_Ns = scale / Ns
    A1_base = zeros(ComplexF64, Ns, (Np + 2))
    build_Abase!(A1_base, s, poles)  # left block
    A1_cplx = zeros(ComplexF64, Ns, Ncols)  # right block
    for n = 1:Nc
        for i = 1:(Np+2)
            A1_cplx[1:Ns, i] .= A1_base[1:Ns, i] .* weight[1:Ns, n]
        end
        for i = 1:Nres
            A1_cplx[1:Ns, (Np+2+i)] .= A1_cplx[1:Ns, i] .* (-f[1:Ns, n])
        end
        A1_reim[1:Ns, :] .= real(A1_cplx)
        A1_reim[(Ns+1):(2Ns), :] .= imag(A1_cplx)
        if relaxed && n == Nc
            A1_reim[end, 1:(Np+2)] .= 0.0
            for i = 1:Nres
                A1_reim[end, Np+2+i] = real(sum(A1_cplx[:, i]) * scale_Ns)
            end
        end
        # Fast VF is a block-wise QR as we only want the last Nres values of the solution. See [3].
        Q, R = qr!(A1_reim)
        i1 = (Np + 3)
        i2 = i1 + Np - 1 + relaxed
        k1 = 1 + (n - 1) * Nres
        k2 = k1 + Np - 1 + relaxed
        A_sys[k1:k2, :] .= R[i1:i2, i1:i2]
        if relaxed && n == Nc
            b_sys[k1:k2] .= Q[end, i1:i2] * scale
        elseif !relaxed
            fw = f[1:Ns, n] .* weight[1:Ns, n]
            b_sys[k1:k2] .= transpose(Q[:, i1:i2]) * [real(fw); imag(fw)]
        end
    end
    norm_cols = [norm(A_sys[:, n]) for n = 1:Nres]
    for n = 1:Nres
        A_sys[:, n] ./= norm_cols[n]
    end
    ldiv!(qr!(A_sys), b_sys)  # x = A \ b
    x_sys[:] ./= norm_cols
    if relaxed
        sig_d = abs(x_sys[end])
        if sig_d < 1e-8
            x_sys[end] = 1e-8 * x_sys[end] / sig_d
            @warn "`d` of sigma too small. Consider stopping execution and setting `relaxed=false`. Resuming..."
        elseif sig_d > 1e8
            x_sys[end] = 1e8 * x_sys[end] / sig_d
            @warn "`d` of sigma too big. Consider stopping execution and setting `relaxed=false`. Resuming..."
        end
        x_sys[1:(end-1)] ./= x_sys[end]  # scale sigma's residues by its `d`
    end
    H = zeros(Np, Np)
    skip_next = false
    for (i, p) in enumerate(poles)
        if skip_next
            skip_next = false
            continue
        elseif isreal(p)
            skip_next = false
            H[:, i] .= -x_sys[i]
            H[i, i] += p
        else
            skip_next = true
            H[1:2:end, i] .= -2.0 * x_sys[i]
            H[i, i] += real(p)
            H[i+1, i] += -imag(p)
            H[1:2:end, i+1] .= -2.0 * x_sys[i+1]
            H[i, i+1] += imag(p)
            H[i+1, i+1] += real(p)
        end
    end
    return eigvals(H)
end


"""
    residue_identification(s, f, poles, weight) -> residues, d, h

Stage 2 of the Vector Fitting. This should be called separately for each column
of `f` and `weight` when `ndims(f) == 2`.

See also [`vector_fitting`](@ref), [`pole_identification`](@ref).
"""
function residue_identification(s, f, poles, weight)
    Ns = length(s)
    Np = length(poles)
    residues = similar(poles)
    Nrows = 2 * Ns
    Ncols = Np + 2
    A1_cplx = zeros(ComplexF64, Ns, Ncols)
    build_Abase!(A1_cplx, s, poles)
    A1_cplx .*= weight
    A_sys = [real(A1_cplx); imag(A1_cplx)]
    norm_cols = [norm(A_sys[:, n]) for n = 1:Ncols]
    for n = 1:Ncols
        A_sys[:, n] ./= norm_cols[n]
    end
    B_cplx = f .* weight
    X_sys = A_sys \ [real(B_cplx); imag(B_cplx)]
    X_sys ./= norm_cols
    skip_next = false
    for (i, p) in enumerate(poles)
        if skip_next
            skip_next = false
            continue
        elseif isreal(p)
            skip_next = false
            residues[i] = X_sys[i]
        else
            skip_next = true
            residues[i] = complex(X_sys[i], X_sys[i+1])
            residues[i+1] = conj(residues[i])
        end
    end
    d = X_sys[Np+1]
    h = X_sys[Np+2]
    return residues, d, h
end


"""
    vector_fitting(
        s,
        f,
        init_poles,
        weight = 1;
        relaxed = true,
        force_stable = true,
        maxiter = 5,
        tol = 1e-12,
    ) -> poles, residues, d, h, fitted, error_norm

Vector Fitting of the array `f` with complex frequency `s`
using a set of initial poles `init_poles`.

`f` can be a matrix of dimensions `(Ns, Nc)` and the fitting will be over
its columns using a set of common poles.

`relaxed` controls the nontriviality constraint. `relaxed=true` usually
converges faster, but can be less stable for non-smooth functions.

`force_stable` controls if unstable poles should be reflected to the semi-left
complex plane, that is, forced to have negative real part.

`maxiter` is the maximum of iterations that will be done to try to achieve a
convergence with desired `error_norm` tolerance `tol`.

See also [`recommended_init_poles`](@ref), [`rational`](@ref),
[`pole_identification`](@ref), [`residue_identification`](@ref).
"""
function vector_fitting(
    s,
    f,
    init_poles,
    weight = 1;
    relaxed = true,
    force_stable = true,
    maxiter = 5,
    tol = 1e-12,
)
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

    if typeof(weight) <: Number
        weight = fill(weight, size(f))
    elseif size(weight) != size(f)
        throw(
            error(
                "It is expected `weight` to be a scalar or to have the same dimensions as f.",
            ),
        )
    end

    poles = sort!(complex(init_poles), by = cplxpair)
    Np = length(poles)
    residues = zeros(ComplexF64, Np, Nc)
    d = zeros(Nc)
    h = zeros(Nc)
    fitted = similar(f)
    local error_norm = Inf
    for iter = 1:maxiter
        if error_norm < tol
            println("convergence achieved at iter. = $(iter)")
            println("error_norm = $(error_norm)")
            break
        end
        poles = pole_identification(s, f, poles, weight, relaxed)
        if force_stable
            for (i, p) in enumerate(poles)
                re_p, im_p = reim(p)
                if re_p > 0.0
                    poles[i] = complex(-re_p, im_p)
                end
            end
        end
        for n = 1:Nc
            residues[:, n], d[n], h[n] =
                residue_identification(s, f[:, n], poles, weight[:, n])
            fitted[:, n] .= rational(s, poles, residues[:, n], d[n], h[n])
        end
        error_norm = norm(f .- fitted, 2)
    end
    perm = sortperm(poles, by = cplxpair)
    poles = poles[perm]
    for n = 1:Nc
        residues[:, n] = residues[perm, n]
    end
    return poles, residues, d, h, fitted, error_norm
end

end  # module
