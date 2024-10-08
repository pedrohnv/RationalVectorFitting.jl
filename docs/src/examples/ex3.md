# [Example 3](@id ex3)

Order 6 fitting of the measured admittance function from a distribution transformer,
which we took from [SINTEF's VFIT3](https://www.sintef.no/en/software/vector-fitting/downloads/vfit3/)

```@example ex3
using RationalVectorFitting
using Plots

fid1 = split.(readlines("03PK10.txt"))
f = zeros(ComplexF64, 160)
for k = 1:160
    A1 = parse(Float64, fid1[k+1][1])
    A2 = parse(Float64, fid1[k+1][2])
    f[k] = A1 * exp(1im * A2 * pi / 180)
end
w = 2pi * range(0, 10e6, length = 401)
w = w[2:161]
s = 1im .* w
freq = imag(s) ./ (2pi)

N = 6  # Order of approximation
init_poles = recommended_init_poles(s, N)

poles, residues, d, h, fitted, error_norm = vector_fitting(s, f, init_poles)

p1 = plot(
    freq,
    abs.(f),
    label = "f(s)",
    linecolor = :blue,
    xlabel = "Frequency [Hz]",
    ylabel = "Magnitude",
    yaxis = :log,
    legend = :right,
)
plot!(freq, abs.(fitted), label = "fitted", linecolor = :darkorange)
plot!(freq, abs.(f - fitted), label = "deviation", linecolor = :green)
savefig(p1, "ex3-mag.svg")

p2 = plot(
    freq,
    rad2deg.(angle.(f)),
    label = "f(s)",
    linecolor = :blue,
    xlabel = "Frequency [Hz]",
    ylabel = "Phase angle [deg]",
    legend = :right,
)
plot!(freq, rad2deg.(angle.(fitted)), label = "fitted", linecolor = :darkorange)
savefig(p2, "ex3-phase.svg")
```

![ex3-mag](ex3-mag.svg)

![ex3-phase](ex3-phase.svg)
