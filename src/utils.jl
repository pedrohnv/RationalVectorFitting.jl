# Utility functions

"""
    cplxpair(x)

To be used to sort an array by real values, then complex conjugate pairs.
The more positive reals will be first, then the pairs with smaller imaginary part.

# Example
```julia
real_values = [-1, -3, -2]
complex_values = [-1 + 1im, -2 - 2im, -1 + 2im]
sorted_values = sort!([complex_values; conj(complex_values); real_values], by = cplxpair)

# output

9-element Vector{Complex{Int64}}:
 -1 + 0im
 -2 + 0im
 -3 + 0im
 -1 - 1im
 -1 + 1im
 -1 - 2im
 -1 + 2im
 -2 - 2im
 -2 + 2im
```
"""
function cplxpair(x)
    return (!isreal(x), abs(imag(x)), abs(real(x)), imag(x))
end
