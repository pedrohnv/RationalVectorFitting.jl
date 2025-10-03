# Elementwise approximation of a 3x3 propagation matrix of an aerial transmission line.
# It corresponds to single propagation mode and time delay is already extracted.

using CSV, DataFrames, LinearAlgebra
using RationalVectorFitting
using Plots

function main()
    # Importing data from a .csv file (lineConstants_H0.csv)
    csv_path = "examples/MODEH_DATA.csv"  # local path! Update with your own path
    # H(w) samples in the frequency domain:
    # Columns organization: index , OMEGA(Ang. frequency), H_00REAL(1st element's real part), H_00IMAG(1st element's imaginary part), ... H_01REAL(2nd element's real part), ...

    # DataFrame with H(w) samples in the frequency domain:
    # Columns organization: index , OMEGA(Ang. frequency), H_00REAL(1st element's real part), H_00IMAG(1st element's imaginary part), ... H_01REAL(2nd element's real part), ...
    Hdata = CSV.read(csv_path, DataFrame)
    w = vec(Hdata[!, "OMEGA"])  # Angular frequency samples
    s = im .* w  # Complex frequency samples
    N = length(w)  # Number of samples

    # Propagation matrix in the frequency domain
    Hw = zeros(ComplexF64, 3, 3, N)

    # Copying data into H:
    k = 3
    for row = 1:3
        for col = 1:3
            real_part = vec(Hdata[!, k])
            imag_part = vec(Hdata[!, k+1])
            Hw[row, col, :] = real_part .+ im .* imag_part  # elements are read in RMO
            k += 2
        end
    end

    # Stacking H(s) data as elements of a frequency domain function F(s).
    # Due to H(s) being asymmetric, F(s) is a flattened version of H(s). Row Major Ordering is used to map H(s) elements into F(s)
    F = zeros(ComplexF64, 3*3, N)
    k = 1  # element index (1-based in Julia)
    for row = 1:3
        for col = 1:3
            F[k, :] = Hw[row, col, :]  # all frequency samples are in z axis
            k += 1
        end
    end
    F = transpose(F)

    n = 35  # Order of approximation

    # Starting poles generation:
    Bet = 10 .^ range(log10(w[1]), log10(w[end]), length = Int(n ÷ 2))
    poles = zeros(ComplexF64, n)

    # setting poles as complex conjugated pairs
    for k = 1:Int(n÷2)
        alf = -Bet[k] / 100
        poles[2*k-1] = alf - im * Bet[k]
        poles[2*k] = alf + im * Bet[k]
    end

    # Using H trace to identify initial poles:
    trH = zeros(ComplexF64, N)
    for k = 1:3
        trH .= trH .+ Hw[k, k, :]
    end

    # for i in 1:10
    #     poles, residues, d, h, fitted, error_norm = vector_fitting(
    #         s,
    #         trH,
    #         poles,
    #         1,
    #         relaxed = false,
    #         force_stable = true,
    #         maxiter = 1,
    #         tol = 1e-12,
    #     )
    # end

    poles, _, _, _, _, _ = vector_fitting(
        s,
        F,
        poles,
        1,
        relaxed = true,
        force_stable = true,
        maxiter = 10,
        tol = 1e-12,
        skip_res = true,
    )

    poles, residues, d, h, fitted, error_norm = vector_fitting(
        s,
        F,
        poles,
        1,
        relaxed = true,
        force_stable = true,
        maxiter = 10,
        tol = 1e-12,
        skip_res = false,
    )

    println("Fitting process completed. Approximation error achieved = ", error_norm)

    begin
        f = w / 2pi
        p1 = plot(
            xaxis = :log,
            ylabel = "Magnitude",
            xlabel = "Frequency [Hz]",
            xticks = exp10.(range(0, 10, 11)),
        )
        for i = 1:size(F, 2)
            if i == 1
                l1 = "F(s)"
                l2 = "fitted(s)"
            else
                l1 = ""
                l2 = ""
            end
            plot!(f, abs.(F[:, i]), label = l1, color = :blue)
            plot!(f, abs.(fitted[:, i]), label = l2, color = :red)
        end
        display(p1)

        p1 = plot(
            xaxis = :log,
            ylabel = "log10(relative error)",
            xlabel = "Frequency [Hz]",
            xticks = exp10.(range(0, 10, 11)),
        )
        for i = 1:size(F, 2)
            rel_error = abs.((F[:, i] - fitted[:, i]) ./ F[:, i])
            e = log10.(rel_error)
            plot!(f, e, label = "", color = :red)
        end
        display(p1)

        max_rel_error = maximum(abs.((F - fitted) ./ F) * 100)
        println("Maximum relative error: $(max_rel_error) %")
    end
end

main()
