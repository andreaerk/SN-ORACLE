using CSV, DataFrames
using NearestNeighbors
using Interpolations
using Statistics
using LinearAlgebra

include("defs.jl")

# Load your data file (replace with your actual path)
fHeS = CSV.read("Data/Models_CC/He_star_summary_Ercolino2025.dat", DataFrame; delim=' ', ignorerepeated=true)

# Extract arrays similar to python
Mc = fHeS.Mconvmax_Heburn
R_20kyr = fHeS.R20kyr
R_10kyr = fHeS.R10kyr
R_5kyr  = fHeS.R5kyr
R_2kyr  = fHeS.R2kyr
R_1kyr  = fHeS.R1kyr
R_1yr   = min.(fHeS.R1yr, 500)
R_max = min.(fHeS.Rmax, 500)

fR_names = [:R20kyr, :R10kyr, :R5kyr, :R2kyr, :R1kyr, :R1yr, :Rmax]
fT = [20e3, 10e3, 5e3, 2e3, 1e3, 1.0, 0.0]  # keep this for z-values
x = Float64[]
y = Float64[]
z = Float64[]

for (i, colname) in enumerate(fR_names)
    for j in 1:length(Mc)
        Rij = fHeS[!, colname][j]
        Rij = min(Rij, 500)  # mimic Python's np.minimum
        push!(x, Mc[j])
        push!(y, log10(Rij))
        push!(z, fT[i])
    end
end

# Build regular grid (similar resolution as python)
xp = range(minimum(x), stop=maximum(x), length=600)
yp = range(minimum(y), stop=maximum(y), length=200)

# Create meshgrid (xp, yp) and interpolate scattered data onto it
# Use Dierckx.jl for scattered interpolation

using Dierckx

# Build spline interpolant for scattered data
spl = Spline2D(x, y, z, kx=1, ky=1, s = 64000000)  # linear spline

# Create zi grid values
zi = [clamp(spl(xv, yv), minimum(z), maximum(z)) for yv in yp, xv in xp]

# Flatten grid points and build KDTree for nearest neighbors in grid and original data
P = hcat(x, y)'
Pi = [ (xp[i], yp[j]) for j in 1:length(yp), i in 1:length(xp)]
Pi_vec = vec(Pi)  # Vector of tuples (Float64, Float64)
Pi_mat = hcat(getindex.(Pi_vec, 1), getindex.(Pi_vec, 2))'  # 2 x N Matrix{Float64}

treeP = KDTree(P)
treePi = KDTree(Pi_mat)

Rmax_null = interpolator(Mc, R_max,    extrapolation=:flat)
Rmax_full = interpolator(Mc, R_20kyr,  extrapolation=:flat)

function retrieve_timepreC(M_convHe, Rl; debug=false)

    if Rl >= Rmax_null(M_convHe) 
        debug && println("too extended -> 0")
        return (0, true)
    elseif Rl <= Rmax_full(M_convHe) 
        debug && println("too small -> 0")
        return (20_000, true)
    end
    p = [M_convHe, log10(Rl)]

    # Query nearest neighbors
    idxP, _ = knn(treeP, [M_convHe; log10(Rl)], 1)    # Query coarse grid
    idxPi, _ = knn(treePi, [M_convHe; log10(Rl)], 1)  # Query fine grid

    zP = z[idxP[1]]    # idxP refers to 'z' (length 77)
    zPi = zi[idxPi[1]] # idxPi refers to 'zi' (length 120000)

    if isnan(zPi)
        # Fallback logic
        if M_convHe > 2.639
            debug && println("overmassive -> 0")
            return 0.0, true
        elseif M_convHe < 1.158
            debug && println("extrapolate from lowest mass")
            return zP, true
        # elseif (M_convHe < 1.6 && Rl < 50) || (M_convHe > 1.6 && Rl < 5)
        #     debug && println("too close -> 20kyr")
        #     return 20_000.0, true
        # elseif (M_convHe < 1.6 && Rl > 100) || (M_convHe > 1.6 && Rl > 4)
        #     debug && println("too extended -> 0")
        #     return 0.0, true
        else
            error("ERROR! CHECK INTERPOLATION!")
        end
    else
        return zPi, false
    end
end

using PythonCall 
using CondaPkg
using Printf
pythondir = "/vol/aibn133/data1/aercolino/SOFTWARE/SN_from_grids/PythonScripts/"
PythonCall.pyimport("sys").path.insert(0, pythondir)
# Import your custom Python script
const interp = pyimport("caseBB_interpolator")
function retrieve_tpreCC_python(m, r)
    try 
        return pyconvert(Float64, interp.retrieve_timepreC(m, r)[0])
    catch 
        return 0
    end
end
print("done\n")

for m in 1:.1:2
    for logr in 0:0.5:3
        r=10^logr
        tP = nothing  
        try 
            tP = retrieve_tpreCC_python(m,r)
        catch 
            continue 
        end
        tN = retrieve_timepreC(m, r)[1]
        rel_diff = (tN-tP)/tP
        @printf("M = %3.1f, R = %5.2e    |||    Python: %5.3e  -  new %5.3e   diff = %.1f%%\n", m, r, tN, tP, rel_diff*100)
    end
end


using CairoMakie
using ColorSchemes


# Define grid
M_vals = 1.:0.05:2.75
R_vals =  (0:0.05:3)
cmap =  reverse(Makie.ColorSchemes.jet)
anchors =  [0,       500,    1e3,    2e3,     5e3,      10e3,   20e3] ./ 20e3
colors =   [:darkgray, :gray, :orange, :red,  :yellow,   :green,  :blue]
# Placeholder: call your Python function here, or load precomputed results
cmap = cgrad(colors, anchors)
# Evaluate on grid
py_results = [(retrieve_tpreCC_python(M, 10^R)) for R in R_vals, M in M_vals]
jl_results = [(retrieve_timepreC(M, 10^R)[1]) for R in R_vals, M in M_vals]  # Assuming returns (value, bool)
difference = (jl_results .- py_results) ./py_results

# Plotting
fig = Figure(resolution=(700, 1000))
ax1 = Axis(fig[1, 1], title="Python results", xlabel="M_convHe", ylabel="Rl")
hm = heatmap!(ax1, M_vals, R_vals, py_results, colormap=cmap,  highclip=:black)
Colorbar(fig[1, 2], hm, label="Time [yr]", scale = log10, )

ax2 = Axis(fig[2, 1], title="Julia results", xlabel="M_convHe", ylabel="Rl", )
hm2 = heatmap!(ax2, M_vals, R_vals, jl_results, colormap=cmap, interpolate=false)
Colorbar(fig[2, 2], hm2, label="Time [yr]")

ax3 = Axis(fig[3, 1], title="Difference (Python - Julia)", xlabel="M_convHe", ylabel="Rl")
hm3 = heatmap!(ax3, M_vals, R_vals, difference, colormap=:balance, colorrange=(-1, 1))
cb = Colorbar(fig[3, 2], hm3, label="Time difference [yr]")
for ax in [ax1,ax2, ax3] 
    lines!(ax, Mc, log10.(R_20kyr), color=:blue)
    lines!(ax, Mc, log10.(R_10kyr), color=:green)
    lines!(ax, Mc, log10.(R_5kyr),  color=:yellow)
    lines!(ax, Mc, log10.(R_2kyr),  color=:brown)
    lines!(ax, Mc, log10.(R_1kyr),  color=:red)
    lines!(ax, Mc, log10.(R_1yr),   color=:gray)
    lines!(ax, Mc, log10.(R_max),   color=:black)
    lines!(ax, Mc, log10.(500) .+ 0 .* Mc,   color=:black)
end

fig

