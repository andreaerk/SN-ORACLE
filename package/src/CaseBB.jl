print("Building Case BB interpolator ...")

# Load your data file (replace with your actual path)
fHeS_BB = CSV.read(data_dir * "Models_CC/He_star_summary_Ercolino2025.dat", DataFrame; delim=' ', ignorerepeated=true)

# Extract arrays similar to python
Mc = fHeS_BB.Mconvmax_Heburn
R_20kyr = fHeS_BB.R20kyr
R_10kyr = fHeS_BB.R10kyr
R_5kyr  = fHeS_BB.R5kyr
R_2kyr  = fHeS_BB.R2kyr
R_1kyr  = fHeS_BB.R1kyr
R_1yr   = min.(fHeS_BB.R1yr, 500)
R_max = min.(fHeS_BB.Rmax, 500)

fR_names = [:R20kyr, :R10kyr, :R5kyr, :R2kyr, :R1kyr, :R1yr, :Rmax]
fT = [20e3, 10e3, 5e3, 2e3, 1e3, 1.0, 0.0]  # keep this for z-values
x = Float64[]
y = Float64[]
z = Float64[]

for (i, colname) in enumerate(fR_names)
    for j in 1:length(Mc)
        Rij = fHeS_BB[!, colname][j]
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

function retrieve_tpreCC(M_convHe, Rl; debug=false)

    if Rl >= Rmax_null(M_convHe) 
        debug && println("too extended -> 0")
        return (0)
    elseif Rl <= Rmax_full(M_convHe) 
        debug && println("too small -> 0")
        return (20_000)
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
            return 0.0
        elseif M_convHe < 1.158
            debug && println("extrapolate from lowest mass")
            return zP
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
        return zPi
    end
end

