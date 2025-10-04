# Current surge fitting in time domain

using LinearAlgebra
using RationalVectorFitting.TimeDomain
using Plots

fid1 = split.(readlines("examples/td_surge.csv"))
N = length(fid1) - 1
t = zeros(N)
vout = zeros(N)
for k = 1:N
    val = split(fid1[k+1][1], ",")
    t[k] = parse(Float64, val[1])
    vout[k] = parse(Float64, val[2])
end

i0 = 31
t0 = t[i0]
t = t[i0:end] .- t0
vout = vout[i0:end]

n = 70  # order
init_poles = -exp10.(range(0, 9, length = n))

dt = t[2] - t[1]
nt = length(t)

# impulse input
vin = zeros(nt)
vin[2] = 1 / dt

niter = 50
has_direct_feedthrough = false
formula = "recursive"
poles = qpol = init_poles
Δt = dt
res = TimeDomain.vector_fitting_time_domain(
    dt,
    vin,
    vout,
    init_poles,
    has_direct_feedthrough = has_direct_feedthrough,
    niter = niter,
    formula = formula,
)
qpol, fitted, pointwise_rmsd, rmsd, pointwise_mean_abs_d, mean_abs_d = res
@show rmsd

begin
    p1 = plot(xlabel = "time [μs]", ylabel = "Magnitude")
    plot!(t*1e6, vout, label = "measured")
    plot!(t*1e6, fitted, label = "fitted")
    display(p1)
end

@show fitted[1]
