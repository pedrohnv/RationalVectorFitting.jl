# [Example 2](@id ex2)

Fitting an order 18 function.

```@example ex2
using RationalVectorFitting
using Plots

poles1 = 2π .* [
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

residues1 = 2π .* [
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

d1 = 0.2
h1 = 2e-5

freq = (range(0, 1e5, length = 200))
s = 2im * pi * freq
f1 = rational(s, poles1, residues1, d1, h1)

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
# now we add the missing half of complex pairs
real_poles = filter(isreal, init_poles)
complex_poles = filter(!isreal, init_poles)
init_poles = sort!([real_poles; complex_poles; conj(complex_poles)], by = cplxpair)

poles, residues, d, h, fitted, error_norm = vector_fitting(s, f1, init_poles)

p1 = plot(
    freq,
    abs.(f1),
    label = "f(s)",
    linecolor = :blue,
    xlabel = "Frequency [Hz]",
    ylabel = "Magnitude",
    yaxis = :log,
    legend = :right,
)
plot!(freq, abs.(fitted), label = "fitted", linecolor = :darkorange)
plot!(freq, abs.(f1 - fitted), label = "deviation", linecolor = :green)
savefig("ex2.svg"); nothing # hide
```

![ex2](ex2.svg)
