#= Time-Domain Vector Fitting

Based on Algorithm 7.3 (page 276) of [1].

[1] Stefano Grivet-Talocia; Bjorn Gustavsen, "The Vector Fitting Algorithm,"
in Passive Macromodeling: Theory and Applications , Wiley, 2016, pp.225-306,
doi: 10.1002/9781119140931.ch7.
=#

module TimeDomainVF

using LinearAlgebra

export convolution, rational_to_state_space, simulate_state_space


@doc raw"""
    rational_to_state_space(poles, residues; real_only=true, reduced=false) -> A, B, C

Converts a pole-residue rational model into a time-domain state-space representation
in canonical Jordan form.

### Arguments
- `poles`: Vector of poles of the rational model.
- `residues`: Vector of residues corresponding to the poles.
- `real_only` (default = `true`):
    If `true`, returns a real-valued state-space system by pairing complex-conjugate poles.
    If `false`, allows the system to be complex.
- `reduced` (default = `false`):
    If `true` reduces the size of the state-space model by:
    1. Discarding half of the conjugate poles.
    2. Doubling the real part of the corresponding residues.
    3. The User discarding the imaginary part during any simulation using this model.

### Returns
- `A`: State matrix.
- `B`: Input matrix.
- `C`: Output matrix.

### Mathematical Representation
The resulting state-space system is represented as:

```math
    \frac{dx(t)}{dt} = A x(t) + B u(t)

    y(t) = C \dot x(t) + D \dot u(t) + E \dot \frac{du(t)}{dt}
```
"""
function rational_to_state_space(poles, residues; real_only = true, reduced = false)
    if ndims(residues) == 1
        nc = 1
        np = length(residues)
    else
        nc, np = size(residues)
    end
    residues = reshape(residues, nc, np)
    if reduced
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
    # TODO passivity constraints
    return am, bm, cm
end


"""
    simulate_state_space(A, B, C, D, input, x0, dt, nt; reduced = false)

Simulation of a state-space model using trapezoidal integration.

### Arguments
- `A`: State matrix of size `(nx, nx)`.
- `B`: Input matrix of size `(nx, n_in)`.
- `C`: Output matrix of size `(n_out, nx)`.
- `D`: Feedthrough of size `(n_out, n_in)`.
- `input`: input matrix of size `(nt, n_in)`.
- `x0`: initial state vector of size `nx`.
- `dt`: time step.
- `nt`: number of time steps.

### Returns
- `output`: output matrix of size `(nt, n_out)`.
"""
function simulate_state_space(A, B, C, D, input, x0, dt, nt)
    # TODO sanity check of the arguments
    Ad = inv(I - dt / 2 * A) * (I + dt / 2 * A)
    Bd = inv(I - dt / 2 * A) * (dt * B)
    if ndims(D) == 0  # D is a scalar
        n_out = n_in = 1
    else
        n_out, n_in = size(D)
    end
    y = zeros(eltype(A), nt, n_out)  # output
    ut = reshape(input, nt, :)  # casts to a column vector if `ndims(input) == 1`
    x = x0  # initialize state vector
    y[1, :] .= C * x + D * ut[1, :]
    for k = 2:nt
        x = Ad * x + Bd * ut[k]  # Update state
        y[k, :] = C * x .+ D * ut[k]  # Compute output
    end
    return y
end


@doc raw"""
    convolution(poles, yt) -> y_j

Calculates the convolution of the time-domain vector `y(t)` (linearly sampled
at time intervals `dt`) with each partial fraction with given `poles`.
Returns the matrix `y_j` of size `(length(yt), length(poles))`.

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
