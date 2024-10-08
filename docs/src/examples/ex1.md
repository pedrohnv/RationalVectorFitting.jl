# [Example 1](@id ex1)

This first example fits an order 3 smooth function using a set of real initial poles.

```@example ex1
using RationalVectorFitting
using Plots

freq = exp10.(range(0, 4, length = 101))  # logspace
s = 2im * pi * freq
poles1 = [-5.0, -100 - 500im, -100 + 500im]
residues1 = [2.0, 30 - 40im, 30 + 40im]
d1 = 0.5
h1 = 0.0
f1 = rational(s, poles1, residues1, d1, h1)

init_poles = -2pi * exp10.(range(0, 4, length = 3))
poles, residues, d, h, fitted, error_norm = vector_fitting(s, f1, init_poles)

p1 = plot(
    freq,
    abs.(f1),
    label = "f(s)",
    linecolor = :blue,
    xlabel = "Frequency [Hz]",
    ylabel = "Magnitude",
    xaxis = :log,
    yaxis = :log,
    legend = :right,
)
plot!(freq, abs.(fitted), label = "fitted", linecolor = :darkorange)
plot!(freq, abs.(f1 - fitted), label = "deviation", linecolor = :green)
savefig("ex1.svg"); nothing # hide
```

![ex1](ex1.svg)
