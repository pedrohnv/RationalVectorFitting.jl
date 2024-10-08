# [Example 4](@id ex4)

Order 50 fitting of the first column of a six-terminal Frequency Domain Network Equivalent (FDNE) system. We took the data from [SINTEF's VFIT3](https://www.sintef.no/en/software/vector-fitting/downloads/vfit3/).

```@example ex4; continued = true
using RationalVectorFitting
using Plots

open("fdne.txt", "r") do fid1
    Nc = parse(Int, readline(fid1))
    Ns = parse(Int, readline(fid1))
    global s = Array{ComplexF64}(undef, Ns)
    global bigY = Array{ComplexF64}(undef, Nc, Nc, Ns)
    for k = 1:Ns
        s[k] = complex(0, parse(Float64, readline(fid1)))
        for row = 1:Nc
            for col = 1:Nc
                a1 = parse(Float64, readline(fid1))
                a2 = parse(Float64, readline(fid1))
                bigY[row, col, k] = complex(a1, a2)
            end
        end
    end
end
f = transpose(bigY[:, 1, :])  # just the 1st column
freq = imag(s) ./ (2pi)

Np = 50  # order of fitting
init_poles = recommended_init_poles(s, Np)
poles, residues, d, h, fitted, error_norm = vector_fitting(s, f, init_poles)
@show error_norm
```

Now we plot the result. Blue lines are the response being fitted, orange lines are the fitted values and green lines are the deviations.

```@example ex4; continued = true
p1 = plot(
    freq,
    abs.(f),
    label = "f(s)",
    linecolor = :blue,
    xlabel = "Frequency [Hz]",
    ylabel = "Magnitude",
    yaxis = :log,
    legend = false,
)
plot!(freq, abs.(fitted), label = "fitted", linecolor = :darkorange)
plot!(freq, abs.(f - fitted), label = "deviation", linecolor = :green)
savefig(p1, "ex4-1.svg")
```

![ex4-1](ex4-1.svg)

That does not look very good. What if we try with a weighting $w(s)$?

```math
w(s) = \frac{1}{\sqrt{f(s)}}
```

```@example ex4; continued = true
weight = @. 1.0 / sqrt(abs(f))
poles, residues, d, h, fitted, error_norm = vector_fitting(s, f, init_poles, weight)
@show error_norm

p1 = plot(
    freq,
    abs.(f),
    label = "f(s)",
    linecolor = :blue,
    xlabel = "Frequency [Hz]",
    ylabel = "Magnitude",
    yaxis = :log,
    legend = false,
)
plot!(freq, abs.(fitted), label = "fitted", linecolor = :darkorange)
plot!(freq, abs.(f - fitted), label = "deviation", linecolor = :green)
savefig(p1, "ex4-2.svg"); nothing # hide
```

![ex4-2](ex4-2.svg)

Damn, it got worse :(

What if we disable the relaxed non-triviality constraint?

```@example ex4
weight = @. 1.0 / sqrt(abs(f))
poles, residues, d, h, fitted, error_norm = vector_fitting(s, f, init_poles, weight; relaxed = false)
@show error_norm

p1 = plot(
    freq,
    abs.(f),
    label = "f(s)",
    linecolor = :blue,
    xlabel = "Frequency [Hz]",
    ylabel = "Magnitude",
    yaxis = :log,
    legend = false,
)
plot!(freq, abs.(fitted), label = "fitted", linecolor = :darkorange)
plot!(freq, abs.(f - fitted), label = "deviation", linecolor = :green)
savefig(p1, "ex4-3.svg")
nothing  # hide
```

![ex4-3](ex4-3.svg)

In some frequency ranges the fitting got better, but worse in others when compared to the first figure. The `error_norm` was smaller.

Unfortunately, the Vector Fitting algorithm relies a lot on trial and error of the user. The fitting could be better if a higher order was used, but that has a tendency to lead to numerically unstable state-space models. Hence, we try to do a low order fitting first.

There are efforts in the literature to make the algorithm more automatic (which some day may be incorporated into this package), see for example:

> A. M. Smith, S. D'Arco, J. A. Suul and B. Gustavsen, "Improved Pole Placement and Compaction of MIMO Vector Fitting Applied to System Identification," in IEEE Transactions on Power Delivery, vol. 39, no. 2, pp. 1259-1270, April 2024, [doi: 10.1109/TPWRD.2024.3364836.](https://doi.org/10.1109/TPWRD.2024.3364836)
