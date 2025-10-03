#= Time-Domain Vector Fitting

Based on Algorithm 7.3 (page 276) of [1].

[1] Stefano Grivet-Talocia; Bjorn Gustavsen, "The Vector Fitting Algorithm,"
in Passive Macromodeling: Theory and Applications , Wiley, 2016, pp.225-306,
doi: 10.1002/9781119140931.ch7.
=#

module TimeDomain

using LinearAlgebra

export convolution,
    rational_to_state_space, symmetric_rational_to_state_space, simulate_state_space


include("utils.jl")

# =============================================================================
# State Space Functions
# =============================================================================

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
```

```math
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

    if length(poles) != np
        throw(error("`poles` must have length `$(np)`"))
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
```

```math
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
    simulate_state_space(A, B, C, D, E, input, dt, nt; x0)

Simulation of an impulsive state-space model using trapezoidal integration.

```math
\frac{dx(t)}{dt} = A x(t) + B u(t)
```

```math
y(t) = C \cdot x(t) + D \cdot u(t) + E \cdot \frac{du(t)}{dt}
```

### Arguments
- `A`: State matrix of size `(nx, nx)`.
- `B`: Input matrix of size `(nx, n_in)`.
- `C`: Output matrix of size `(n_out, nx)`.
- `D`: Feedthrough of size `(n_out, n_in)`.
- `E`: Impulsive term of size `(n_out, n_in)`.
- `input`: input matrix `u(t)` of size `(nt, n_in)`.
- `dt`: time step.
- `nt`: number of time steps including.
- `x0` (optional): initial state vector of size `(nx,)`. Default is zero.

### Returns
- `output`: output matrix `y(t + dt/2)` of size `(nt, n_out)`.
    Note that `output[1]` corresponds to `y(t = dt/2; x = x0)` because of the implicit scheme.
"""
function simulate_state_space(A, B, C, D, E, input, dt, nt, x0 = nothing)
    if ndims(A) == 1
        nx = size(A, 1)
        A = reshape(A, nx, nx)
    end

    if ndims(B) == 1
        B = reshape(B, length(B), 1)
    end

    if ndims(C) == 1
        C = reshape(C, 1, length(C))
    end

    if ndims(D) == 1
        D = reshape(D, 1, length(D))
    end

    if ndims(E) == 1
        E = reshape(E, 1, length(E))
    end

    nx = size(B, 1)
    n_out, n_in = size(D)

    if x0 === nothing
        x0 = zeros(nx)
    end

    if length(x0) < nx
        throw(error("Initial state vector must have `length(x0) ≥ $(nx)`"))
    end

    inv_I_A_dt_2 = inv(I - A * dt / 2)
    Ar = inv_I_A_dt_2 * (I + A * dt / 2)
    Br0 = inv_I_A_dt_2 * B * (dt / 2)
    Br1 = Br0  # because trapezoidal integration
    # Do a state-variable transformation `xr := x + Br0 * u`
    Br1 = Br1 + Ar * Br0
    Dr = D .+ C * Br0
    Cr = C

    if maximum(abs.(E)) > 0
        Ar = cat(Ar, -I(n_in), dims = (1, 2))
        Cr = [Cr I(n_in)]
        Br1 = [Br1; -4 * E / dt]
        Dr .+= 2 * E / dt
        x = zeros(ComplexF64, nx + n_in)
        # Set initial condition for extended state
        if length(x0) == nx
            x[1:nx] = x0
        else
            x = x0
        end
    else
        x = x0
    end
    u = reshape(input, nt, :)  # casts to a column vector if `ndims(input) == 1`
    y = zeros(eltype(A), nt, n_out)
    # `k = 1` corresponds to `t = 0`
    for k = 1:nt
        y[k, :] .= Cr * x + Dr .* u[k, :]
        x[:] .= Ar * x + Br1 .* u[k, :]  # in the next step
    end
    return y
end


@doc raw"""
    convolution_terms(dt, yt, poles, formula = "recursive") -> y_j

Calculates the convolution of the time-domain vector `y(t)` (linearly sampled
at time intervals `dt`) with each partial fraction with given `poles`.
Returns the matrix `y_j` of size `(length(yt), length(poles))`.

The formula used can be selected with the `formula` argument. Options are
"recursive" (default) or "trapezoidal".

The convolution result can be obtained with
    `[sum(y_j[k, :] .* residues) for k = 1:nt]`

```math
\sum_{j = 1}^{N} \mathcal{L}^{-1} \left\{ \frac{1}{s - a_n} \right\} * y(t)
```

See also [`convolution`](@ref).
"""
function convolution_terms(dt, yt, poles, formula = "recursive")
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
    return y_j
end


@doc raw"""
    convolution(dt, yt, poles, residues, formula = "recursive") -> x

Calculates the convolution of the time-domain vector `y(t)` (linearly sampled
at time intervals `dt`) with each partial fraction with given `poles` and
`residues`. Returns the convoluted signal:

```math
x(t) = \sum_{j = 1}^{N} \mathcal{L}^{-1} \left\{ \frac{1}{s - a_n} \right\} * y(t)
```

The formula used can be selected with the `formula` argument. Options are
"recursive" (default) or "trapezoidal".

See also [`convolution_terms`](@ref).
"""
function convolution(dt, yt, poles, residues, formula = "recursive")
    y_j = convolution_terms(dt, yt, poles, formula)
    return [sum(y_j[k, :] .* residues) for k = 1:nt]
end


# =============================================================================
# Vector Fitting
# =============================================================================


"""Identify complex conjugate pairs."""
function idxLine(polos)
    npol = length(polos)
    indices = zeros(Int, npol)
    for mm = 1:npol
        if !isreal(polos[mm])
            if mm == 1
                indices[mm] = 1
            elseif indices[mm-1] == 0 || indices[mm-1] == 2
                indices[mm] = 1
                if mm + 1 <= npol
                    indices[mm+1] = 2
                end
            else
                indices[mm] = 2
            end
        end
    end
    return indices
end


"""Process complex lines for real-valued matrix output."""
function zAuxCompLine(z1, index)
    npol = length(z1)
    result = zeros(npol)
    for mm = 1:npol
        if index[mm] == 0
            result[mm] = real(z1[mm])
        elseif index[mm] == 1
            result[mm] = 2 * real(z1[mm])
        elseif index[mm] == 2
            result[mm] = 2 * imag(z1[mm])
        end
    end
    return result
end


"""Calculate coefficient matrix A for the system."""
function coeffAA(
    Δt,
    vin,
    vout,
    poles;
    has_direct_feedthrough = false,
    formula = "recursive",
)
    npol = length(poles)
    z1 = convolution_terms(Δt, vin, poles, formula)
    z2 = convolution_terms(Δt, vout, poles, formula)

    # Remove small values
    z1 = [abs(z) < 1e-10 ? 0.0 : z for z in z1]
    z2 = [abs(z) < 1e-10 ? 0.0 : z for z in z2]

    nt = length(vin)
    AAAux = zeros(ComplexF64, nt-1, 2*npol + has_direct_feedthrough)
    AAtemp = zeros(nt-1, 2*npol + has_direct_feedthrough)

    if has_direct_feedthrough
        AAAux[:, 1:npol] = z1[2:end, :]
        AAAux[:, npol+1] = vin[2:end]
        AAAux[:, (npol+2):(2*npol+1)] = -z2[2:end, :]
    else
        AAAux[:, 1:npol] = z1[2:end, :]
        AAAux[:, (npol+1):(2*npol)] = -z2[2:end, :]
    end

    # Process complex lines to real-valued matrix
    for nm = 1:(nt-1)
        indicesz1 = idxLine(AAAux[nm, :])
        AAtemp[nm, :] = zAuxCompLine(AAAux[nm, :], indicesz1)
    end

    return AAtemp
end


"""Calculate coefficient matrix T for the system."""
function coeffTT(Δt, vin, poles; has_direct_feedthrough = false, formula = "recursive")
    npol = length(poles)
    z1 = convolution_terms(Δt, vin, poles, formula)

    # Remove small values
    z1 = [abs(z) < 1e-10 ? 0.0 : z for z in z1]

    nt = length(vin)
    AAAux = zeros(ComplexF64, nt-1, npol + has_direct_feedthrough)
    AAtemp = zeros(nt-1, npol + has_direct_feedthrough)

    if has_direct_feedthrough
        AAAux[:, 1:npol] = z1[2:end, :]
        AAAux[:, npol+1] = vin[2:end]
    else
        AAAux[:, 1:npol] = z1[2:end, :]
    end

    # Process complex lines to real-valued matrix
    for nm = 1:(nt-1)
        indicesz1 = idxLine(AAAux[nm, :])
        AAtemp[nm, :] = zAuxCompLine(AAAux[nm, :], indicesz1)
    end

    return AAtemp
end


"""Time-Domain Vector Fitting.

### Arguments
- Δt: time step [s]
- vin: input vector `x(t)`. It is assumed that `vin[1] = 0`.
- vout: output vector `y(t)`. It is assumed that `vout[1] = 0`.
- init_poles: initial poles guess of `H(s)`
- has_direct_feedthrough: bool to include or not a direct feedthrough term
- niter: number of iterations
- formula: to use in the convolution, either "trapezoidal" or "recursive"

### Returns
- new poles
- fitted values of `y_fit(t)`
- pointwise root mean squared difference
- root mean squared difference
- pointwise mean absolute difference
- mean absolute difference
"""
function vector_fitting_time_domain(
    Δt,
    vin,
    vout,
    init_poles;
    has_direct_feedthrough = false,
    niter = 5,
    formula = "recursive",
)
    if niter < 1
        throw(ArgumentError("niter must be greater or equal to 1"))
    end

    npol = length(init_poles)
    nt = length(vin)

    qpol = copy(init_poles)
    t = range(0, length = nt) * Δt

    local qpol, fitted, pointwise_rmsd, rmsd, pointwise_mean_abs_d, mean_abs_d
    for nn = 1:niter
        # Calculate coefficient matrix A
        AAtemp = coeffAA(
            Δt,
            vin,
            vout,
            qpol,
            has_direct_feedthrough = has_direct_feedthrough,
            formula = formula,
        )
        BBtemp = vout[2:end]

        # Get indices for complex conjugate pairs
        indices = idxLine(qpol)

        # QR decomposition
        Q, R = qr(AAtemp)
        m, n = size(AAtemp)
        k = min(m, n)
        Q = Q[:, 1:k]
        R = R[1:k, 1:k]
        AAA = R[(end-npol+1):end, (end-npol+1):end]
        BBB = transpose(Q[:, (end-npol+1):end]) * BBtemp

        # Scale the matrix
        normas = norm.(eachcol(AAA), 2)
        AAAaux = AAA ./ transpose(normas)
        Xaux = AAAaux \ BBB
        scale = 1.0 ./ normas
        XX = real(Xaux .* scale)

        # Build system matrices
        Azeros = zeros(npol, npol)
        Bzeros = zeros(npol, 1)
        Czeros = XX[1:npol]

        # Set up Bzeros based on indices
        for mm = 1:npol
            if indices[mm] == 0
                Bzeros[mm] = 1.0
            elseif indices[mm] == 1
                Bzeros[mm] = 2.0
            else
                Bzeros[mm] = 0.0
            end
        end

        # Build Azeros matrix
        for mm = 1:npol
            if indices[mm] == 0
                Azeros[mm, mm] = real(qpol[mm])
            elseif indices[mm] == 1
                Azeros[mm, mm] = real(qpol[mm])
                Azeros[mm, mm+1] = imag(qpol[mm])
                Azeros[mm+1, mm] = -imag(qpol[mm])
                Azeros[mm+1, mm+1] = real(qpol[mm])
            end
        end

        # Calculate eigenvalues
        polaux = eigvals(Azeros - Bzeros * transpose(Czeros))
        polaux = [abs(p) < 1e-10 ? 0.0 : p for p in polaux]  # Equivalent to Chop
        qpol = -abs.(real.(polaux)) .+ 1im .* imag.(polaux)

        # Calculate residues
        indicesffT = idxLine(qpol)

        AAtemp = coeffTT(
            Δt,
            vin,
            qpol,
            has_direct_feedthrough = has_direct_feedthrough,
            formula = formula,
        )
        BBtemp = vout[2:end]

        resfim = zeros(ComplexF64, npol)

        # QR decomposition for residue calculation
        Q, R = qr(AAtemp)
        m, n = size(AAtemp)
        k = min(m, n)
        Q = Q[:, 1:k]
        R = R[1:k, 1:k]
        AAA = R[
            (end-npol-has_direct_feedthrough+1):end,
            (end-npol-has_direct_feedthrough+1):end,
        ]
        BBB = transpose(Q[:, (end-npol-has_direct_feedthrough+1):end]) * BBtemp

        # Scale the matrix
        normas = norm.(eachcol(AAA))
        AAAaux = AAA ./ transpose(normas)
        Xaux = AAAaux \ BBB
        scale = 1.0 ./ normas
        XXff = real(Xaux .* scale)

        # Calculate residues
        i = 1
        while i <= npol
            if indicesffT[i] == 0
                resfim[i] = XXff[i]
                i += 1
            elseif indicesffT[i] == 1
                resfim[i] = XXff[i] + 1im * XXff[i+1]
                resfim[i+1] = XXff[i] - 1im * XXff[i+1]
                i += 2
            else
                i += 1
            end
        end

        if has_direct_feedthrough
            dd = XXff[end]
        else
            dd = 0.0
        end

        # Calculate output function
        function yout(t_val)
            result = dd
            for n = 1:npol
                result += resfim[n] * exp(t_val * qpol[n])
            end
            return result
        end

        fitted = real.(yout.(t))

        # Calculate error metrics
        pointwise_rmsd = sqrt.(abs2.(vout - fitted) ./ nt)
        rmsd = sqrt(sum(abs2.(vout - fitted)) / nt)
        pointwise_mean_abs_d = abs.(vout - fitted) ./ nt
        mean_abs_d = sum(abs.(vout - fitted)) / nt
    end

    return qpol, fitted, pointwise_rmsd, rmsd, pointwise_mean_abs_d, mean_abs_d
end

end  # Module
