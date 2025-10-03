# Current surge fitting in time domain

using LinearAlgebra
using RationalVectorFitting.TimeDomain

fid1 = split.(readlines("examples/td_surge.csv"))
N = length(fid1) - 1
t = zeros(N)
vout = zeros(N)
for k = 1:N
    val = split(fid1[k+1][1], ",")
    t[k] = parse(Float64, val[1])
    vout[k] = parse(Float64, val[2])
end

t0 = t[1]
t = t[31:end] .- t0
vout = vout[31:end]

n = 70  # order
init_poles = -exp10.(range(0, 9, length = n))

dt = t[2] - t[1]
nt = length(t)

# impulse input
vin = zeros(nt)
vin[2] = 1 / dt
#vin .+= 1e-12  # to avoid 1/0

niter = 50
has_direct_feedthrough = true
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

using Plots
begin
    p1 = plot(t, vout)
    plot!(t, fitted)
    display(p1)
end
