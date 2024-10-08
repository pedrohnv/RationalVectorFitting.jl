# RationalVectorFitting

Fast Relaxed Vector Fitting implementation in Julia.

Given a transfer function $f(s) = y$, the Vector Fitting algorithm tries to find a rational approximation

$$f(s) \approx \sum_{n=1}^N \frac{r_n}{s - a_n} + d + s h$$

where $s$ is the complex frequency, $r_n$ are the complex residues, $a_n$ are the complex poles, $d$ and $h$ are real constants.

The transfer function can be a vector $f(s) = \[y_1, \dots, y_m\]$ and the Vector Fitting algorithm will fit the response using the same set of poles $a_n$ for all $y_m$.

A rational representation of a transfer function makes it easier to find a [state space canonical realization](https://en.wikipedia.org/wiki/Realization_(systems)#Canonical_realizations) of a system and to [perform convolutions](https://doi.org/10.4236/jamp.2022.106144).

## Usage Example

```julia
using RationalVectorFitting
using Plots

Ns = 101
freq = exp10.(range(0, 4, length = Ns))
s = 2im * pi * freq
poles0 = [-5.0, -100 - 500im, -100 + 500im]
residues0 = [2.0, 30 - 40im, 30 + 40im]
d0 = 0.5
h0 = 0.0
f = rational(s, poles0, residues0, d0, h0)
init_poles = -2pi * exp10.(range(0, 4, length = 3))
poles, residues, d, h, fitted, error_norm = vector_fitting(s, f, init_poles)
begin
    p1 = plot(freq, abs.(f), label="f(s)", linecolor=:blue, xlabel="Frequency [Hz]", xaxis=:log, yaxis=:log, legend=:right)
    plot!(freq, abs.(fitted), label="fitted(s)", linecolor=:darkorange)
    plot!(freq, abs.(f - fitted), label="deviation", linecolor=:green)
    display(p1)
end
```

---

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://pedrohnv.github.io/RationalVectorFitting.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://pedrohnv.github.io/RationalVectorFitting.jl/dev)
[![Build Status](https://github.com/pedrohnv/RationalVectorFitting.jl/workflows/Test/badge.svg)](https://github.com/pedrohnv/RationalVectorFitting.jl/actions)
[![Test workflow status](https://github.com/pedrohnv/RationalVectorFitting.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/pedrohnv/RationalVectorFitting.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/pedrohnv/RationalVectorFitting.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/pedrohnv/RationalVectorFitting.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/pedrohnv/RationalVectorFitting.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/pedrohnv/RationalVectorFitting.jl/actions/workflows/Docs.yml?query=branch%3Amain)

[![Coverage](https://codecov.io/gh/pedrohnv/RationalVectorFitting.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/pedrohnv/RationalVectorFitting.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/pedrohnv/RationalVectorFitting.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

## How to Cite

If you use RationalVectorFitting.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/pedrohnv/RationalVectorFitting.jl/blob/main/CITATION.cff).

## Contributing

If you want to make contributions of any kind, please first that a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://pedrohnv.github.io/RationalVectorFitting.jl/dev/90-contributing/).

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

## Bibliography

[1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses by vector fitting," in IEEE Transactions on Power Delivery, vol. 14, no. 3, pp. 1052-1061, July 1999, [doi: 10.1109/61.772353](https://doi.org/10.1109/61.772353).

[2] B. Gustavsen, "Improving the pole relocating properties of vector fitting," in IEEE Transactions on Power Delivery, vol. 21, no. 3, pp. 1587-1592, July 2006, [doi: 10.1109/TPWRD.2005.860281](https://doi.org/10.1109/TPWRD.2005.860281).

[3] D. Deschrijver, M. Mrozowski, T. Dhaene and D. De Zutter, "Macromodeling of Multiport Systems Using a Fast Implementation of the Vector Fitting Method," in IEEE Microwave and Wireless Components Letters, vol. 18, no. 6, pp. 383-385, June 2008, [doi: 10.1109/LMWC.2008.922585](https://doi.org/10.1109/LMWC.2008.922585).

[4] A. M. Smith, S. D'Arco, J. A. Suul and B. Gustavsen, "Improved Pole Placement and Compaction of MIMO Vector Fitting Applied to System Identification," in IEEE Transactions on Power Delivery, vol. 39, no. 2, pp. 1259-1270, April 2024, [doi: 10.1109/TPWRD.2024.3364836.](https://doi.org/10.1109/TPWRD.2024.3364836)
