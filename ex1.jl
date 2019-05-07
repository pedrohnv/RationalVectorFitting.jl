# ex1.m
#
# Fitting an artificially created frequency response (single element)
#
# -Creating a 3rd order frequency response f(s)
# -Fitting f(s) using vectfit3.m
#   -Initial poles: 3 logarithmically spaced real poles
#   -1 iteration
#
# This example script is part of the vector fitting package (VFIT3.zip)
# Last revised: 08.08.2008.
# Created by:   Bjorn Gustavsen.
# ================================
#include("vectfit3.jl")
# logspace = exp10.(range(start, stop, length))
#Frequency samples:
Ns = 101;
s = 2*pi*1im*exp10.(range(0, stop=4, length=Ns));

f = zeros(ComplexF64, 1, Ns);
for k=1:Ns
    global sk = s[k];
    global f[1,k] = 2/(sk+5) + (30+40im)/(sk-(-100+500im)) + (30-40im)/(sk-(-100-500im)) + 0.5;
end

#Initial poles for Vector Fitting:
N = 3; #order of approximation
poles = Array{ComplexF64}( -2*pi*exp10.(range(0, stop=4, length=N)) ); #Initial poles

weight=ones(1,Ns);
relax=true;
stable=true;
asymp=3;
skip_pole=false;
skip_res=false;

poles, rmserr, fit = vectfit3(f, s, poles, weight, relax, stable, asymp,
                              skip_pole, skip_res)
f = transpose(f);
fit = transpose(fit);
println("rmserr = ", rmserr)
using Plots
plotly()
# Abs
plot(imag(s), map(abs, f), xaxis=:log, yaxis=:log, label="func.",
     ylabel="Magnitude [p.u.]", xlabel="freq. (Hz)")
plot!(imag(s), map(abs, fit), label="fitted")
plot!(imag(s), map(abs, f - fit), label="error")
# Phase
plot(imag(s), map(rad2deg, map(angle, f)), xaxis=:log, label="func.",
     ylabel="Phase [deg]", xlabel="freq. (Hz)")
plot!(imag(s), map(rad2deg, map(angle, fit)), label="fitted")
#plot!(imag(s), map(rad2deg, map(angle, f - fit)), label="error")
