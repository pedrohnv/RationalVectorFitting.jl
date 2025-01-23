#= Time-Domain Vector Fitting

Based on Algorithm 7.3 (page 276) of [1].

[1] Stefano Grivet-Talocia; Bjorn Gustavsen, "The Vector Fitting Algorithm,"
in Passive Macromodeling: Theory and Applications , Wiley, 2016, pp.225-306,
doi: 10.1002/9781119140931.ch7.
=#

module TimeDomainVF

using LinearAlgebra

export convolution,
    rational_to_state_space, symmetric_rational_to_state_space, simulate_state_space


@doc raw"""
    rational_to_state_space(poles, residues; real_only=true, reduced=false) -> A, B, C

Converts a pole-residue rational model into a time-domain impulsive state-space
representation in canonical Jordan form.

### Arguments
- `poles`: Vector containing the `np` set of common poles of the rational model.
- `residues`: Vector (size `(np,)`) or Matrix (size `(nc, np)`) of residues corresponding to the poles.
- `real_only` (default = `true`):
    If `true`, the function returns a real-valued state-space system by pairing complex-conjugate poles.
    If `false`, the system may contain complex-valued matrices.
- `reduced` (default = `false`):
    If `true`, the function returns a reduced complex state-space model by:
    1. Discarding half of the conjugate poles.
    2. Doubling the real part of the corresponding residues.
    3. Users should discard the imaginary part during simulations when using this model.
    **Note:** This option overrides the `real_only` parameter.

### Returns
- `A`: State matrix.
- `B`: Input matrix.
- `C`: Output matrix.

### Mathematical Representation
The resulting impulsive state-space system is given by:

```math
    \frac{dx(t)}{dt} = A x(t) + B u(t)

    y(t) = C \dot x(t) + D \dot u(t) + E \dot \frac{du(t)}{dt}
```

See also [`symmetric_rational_to_state_space`](@ref).
"""
function rational_to_state_space(poles, residues; real_only = true, reduced = false)
    # TODO maybe it would be better to use Sparse Arrays?
    # TODO enforce passivity constraints?
    if ndims(residues) == 1
        nc = 1
        np = length(residues)
    elseif ndims(residues) == 2
        nc, np = size(residues)
    else
        throw(error("It was expected that `residues` had ndims equal to 1 or 2."))
    end
    residues = reshape(residues, nc, np)
    if reduced
        if real_only
            @warn "`real_only` was ignored because `reduced` is true"
        end
        idx = imag(poles) .>= 0.0
        halfpoles = poles[idx]
        np = length(halfpoles)
        am = diagm(halfpoles)
        bm = ones(np)
        cm = residues[:, idx]
        for i = 1:np
            if !isreal(halfpoles[i])
                cm[:, i] .*= 2.0
            end
        end
        return am, bm, cm
    end

    if real_only
        am = zeros(np, np)
        bm = zeros(np)
        cm = zeros(nc, np)
        i = 1
        while i <= np
            if isreal(poles[i])
                am[i, i] = real(poles[i])
                bm[i] = 1.0
                cm[:, i] .= real(residues[:, i])
                i += 1
            else
                am[i, i] = real(poles[i])
                am[i+1, i] = -imag(poles[i])
                am[i, i+1] = imag(poles[i])
                am[i+1, i+1] = real(poles[i])
                bm[i] = 2.0
                bm[i+1] = 0.0
                cm[:, i] .= real(residues[:, i])
                cm[:, i+1] .= imag(residues[:, i])
                i += 2
            end
        end
    else
        am = zeros(ComplexF64, np, np)
        bm = zeros(ComplexF64, np)
        cm = zeros(ComplexF64, nc, np)
        for i = 1:np
            am[i, i] = poles[i]
            bm[i] = 1.0
            cm[:, i] .= residues[:, i]
        end
    end
    return am, bm, cm
end


@doc raw"""
    symmetric_rational_to_state_space(poles, residues; real_only=true, reduced=false) -> A, B, C

Like [`rational_to_state_space`](@ref), converts a pole-residue rational model
into a time-domain impulsive state-space representation in canonical Jordan form.

### Arguments
- `poles`: Vector containing the `np` set of common poles of the rational model.
- `residues`: Symmetric matrix array of size `(nc, nc, np)` of the residues corresponding to the poles.
- `real_only` (default = `true`):
    If `true`, the function returns a real-valued state-space system by pairing complex-conjugate poles.
    If `false`, the system may contain complex-valued matrices.
- `reduced` (default = `false`):
    If `true`, the function returns a reduced complex state-space model by:
    1. Discarding half of the conjugate poles.
    2. Doubling the real part of the corresponding residues.
    3. Users should discard the imaginary part during simulations when using this model.
    **Note:** This option overrides the `real_only` parameter.

### Returns
- `A`: State matrix.
- `B`: Input matrix.
- `C`: Output matrix.

### Mathematical Representation
The resulting impulsive state-space system is given by:

```math
    \frac{dx(t)}{dt} = A x(t) + B u(t)

    y(t) = C \cdot x(t) + D \cdot u(t) + E \cdot \frac{du(t)}{dt}
```

See also [`rational_to_state_space`](@ref).
"""
function symmetric_rational_to_state_space(
    poles,
    residues;
    real_only = true,
    reduced = false,
)
    npoles = length(poles)
    nc1, nc2, npr = size(residues)
    if npoles != npr
        throw(error("`length(poles) != size(residues)[3]`"))
    end
    if nc1 != nc2 || !all([issymmetric(residues[:, :, i]) for i = 1:npoles])
        throw(error("`residues[:, :, i]` is not symmetric"))
    end
    # flatten half of the system
    nh = Int(nc1 * (nc1 + 1) / 2)
    nc = Int((sqrt(1 + 8 * nh) - 1) / 2)
    res = zeros(ComplexF64, (nh, npoles))
    let nr = 1
        for i = 1:nc
            for k = i:nc
                res[nr, :] .= residues[k, i, :]
                nr += 1
            end
        end
    end
    av, bv, cv =
        rational_to_state_space(poles, res; real_only = real_only, reduced = reduced)
    # Recover the symmetric matrix characteristic of the system
    np = size(cv)[2]
    am = cat([av for i = 1:nc]..., dims = (1, 2))
    bm = zeros(nc * np, nc)
    cm = zeros(ComplexF64, nc, nc * np)
    let nr = 1
        for k = 1:nc
            k1 = 1 + np * (k - 1)
            k2 = k1 + np - 1
            bm[k1:k2, k] .= bv
            for i = k:nc
                i1 = 1 + np * (i - 1)
                i2 = i1 + np - 1
                cm[i, k1:k2] .= cv[nr, :]
                cm[k, i1:i2] .= cv[nr, :]
                nr += 1
            end
        end
    end
    return am, bm, cm
end


@doc raw"""
    simulate_state_space(A, B, C, D, E, input, x0, dt, nt)

Simulation of an impulsive state-space model using trapezoidal integration.

```math
    \frac{dx(t)}{dt} = A x(t) + B u(t)

    y(t) = C \cdot x(t) + D \cdot u(t) + E \cdot \frac{du(t)}{dt}
```

### Arguments
- `A`: State matrix of size `(nx, nx)`.
- `B`: Input matrix of size `(nx, n_in)`.
- `C`: Output matrix of size `(n_out, nx)`.
- `D`: Feedthrough of size `(n_out, n_in)`.
- `E`: Impulsive term of size `(n_out, n_in)`.
- `input`: input matrix of size `(nt, n_in)`.
- `dt`: time step.
- `nt`: number of time steps.

### Returns
- `output`: output matrix of size `(nt, n_out)`.
"""
function simulate_state_space(A, B, C, D, E, input, dt, nt)
    # TODO sanity check of the arguments
    inv_I_A_dt_2 = inv(I - A * dt / 2)
    Ar = inv_I_A_dt_2 * (I + A * dt / 2)
    Br0 = inv_I_A_dt_2 * B * (dt / 2)
    Br1 = Br0  # because trapezoidal integration
    # Do a state-variable transformation `xr := x + Br0 * u`
    Br1 = Br1 + Ar * Br0
    Dr = D .+ C * Br0
    Cr = C
    nx = size(B)[1]
    if ndims(D) == 0  # D is scalar
        n_in = n_out = 1
    elseif ndims(D) == 1
        n_out = size(D)[1]
        n_in = 1
    else
        n_out, n_in = size(D)
    end
    if maximum(abs.(E)) > 0
        Ar = cat(Ar, -I(n_in), dims = (1, 2))
        Cr = [Cr I(n_in)]
        Br1 = [Br1; -4 * E / dt]
        Dr .+= 2 * E / dt
        x = zeros(ComplexF64, nx + n_in)
    else
        x = zeros(ComplexF64, nx)
    end
    u = reshape(input, nt, :)  # casts to a column vector if `ndims(input) == 1`
    y = zeros(eltype(A), nt, n_out)
    for k = 1:nt
        y[k, :] .= Cr * x + Dr .* u[k, :]
        x[:] .= Ar * x + Br1 .* u[k, :]  # in the next step
    end
    return y
end


@doc raw"""
    convolution(poles, residues, yt, dt, formula = "recursive") -> y_j

Calculates the convolution of the time-domain vector `y(t)` (linearly sampled
at time intervals `dt`) with each partial fraction with given `poles` and
`residues`. Returns the matrix `y_j` of size `(length(yt), length(poles))`.

The formula used can be selected with the `formula` argument. Options are
"recursive" (default) or "trapezoidal".

```math
\sum_{j = 1}^{N} \mathcal{L}^{-1} \left\{ \frac{1}{s - a_n} \right\} * y(t)
```
"""
function convolution(poles, residues, yt, dt, formula = "recursive")
    if formula == "recursive"
        qn_dt = poles .* dt
        qn2_dt = poles .* qn_dt
        alphaj = exp.(qn_dt)
        betaj0 = @. (-1 - qn_dt + alphaj) / qn2_dt
        betaj1 = @. (1 + (qn_dt - 1) * alphaj) / qn2_dt
    elseif formula == "trapezoidal"
        qn_dt_2 = poles .* dt / 2
        den = (1 .- qn_dt_2)
        alphaj = @. (1 + qn_dt_2) / den
        betaj0 = @. (dt / 2) / den
        betaj1 = betaj0
    else
        throw(
            error(
                """Unknown option for `formula`. It must be "recursive" or "trapezoidal".""",
            ),
        )
    end
    nt = length(yt)
    np = length(poles)
    y_j = zeros(ComplexF64, nt, np)
    for i = 1:np
        for k = 2:nt
            y_j[k, i] = alphaj[i] * y_j[k-1, i] + betaj0[i] * yt[k] + betaj1[i] * yt[k-1]
        end
    end
    return [sum(y_j[k, :] .* residues) for k = 1:nt]
end


# TODO time-domain vector fitting

end  # module
