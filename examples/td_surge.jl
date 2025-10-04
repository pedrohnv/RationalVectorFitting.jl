# Current surge fitting in time domain

using LinearAlgebra
using RationalVectorFitting.TimeDomain

fid1 = split.(readlines("examples/surto1_filtered.csv"))
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
using LaTeXStrings

# Set up LaTeX-style fonts and publication quality settings
   default(
       fontfamily = "Computer Modern",
       titlefontsize = 14,
       guidefontsize = 16,
       tickfontsize = 10,
       legendfontsize = 12,
       linewidth = 2,
       framestyle = :box,
       grid = true,
       gridlinewidth = 0.5,
       gridalpha = 0.3,
       minorgrid = true,
       minorgridlinewidth = 0.25,
       minorgridalpha = 0.15,
       dpi = 600,
       margin = 5Plots.mm,
       size = (800, 450),
    )
# plot results
gr()

begin
    # Define colors for better visibility
     color_data = :black
     color_fitted = :orange
     color_error = :cyan
    p1 = plot(t * 1e6, vout,
         color = color_data,
         legend = :topright,
         linewidth = 2,
         linestyle = :solid,
         label = "data",  
         #title = "Time Domain Vector Fitting"
         alpha = 0.8)
    plot!(t * 1e6, fitted,
         color=color_fitted,
         linewidth = 2,
         linestyle = :dash, 
         legend = :topright, 
         label = "fitted",
         alpha = 0.9)
    plot!(t * 1e6, vout - fitted, 
          color = color_error,
          linewidth = 2,
          linestyle = :dot,
          alpha = 0.7,
          legend = :topright, label = "error")

    # Set labels with LaTeX formatting
      xlabel!(L"\mathrm{Time}~(\mu\mathrm{s})")
      ylabel!(L"\mathrm{Current}~(\mathrm{kA})")
      title!(L"\mathrm{Time~Domain~Vector~Fitting}")
  
    # Configure legend
      plot!(
           legend = :topright,
           legendfontsize = 10,
           legendtitle = nothing,
           background_color_legend = :white,
           foreground_color_legend = :black,
           legend_font_color = :black
     )
 
    # Set size for publication
    #  plot!(size = (1600, 900))
  
    # Display the plot
    display(p1)
end


