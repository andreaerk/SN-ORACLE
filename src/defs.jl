## defs.jl


G = 6.6743e-8 #g-1 cm3 s-1
L☼ = 3.839e33 #egs s-1
M☼ = 1.989e33 #g
R☼ = 6.955e10 #cm
yr = 3.154e+7 #s
km   = 1e5 #cm
kms  = 1e5 #cm/s
day = 24*60*60 #s
kB = 1.380649e-16 #erg/s
amu = 1.6605402e-24 #g
sigma_boltz = 5.6704e-5 #gs-3K-4
c = 2.99792458e10 #cm/s
c_light = 2.99792458e10 #cm/s
a = 4*sigma_boltz/c # gs⁻³K⁻⁴ / (cm/s) = g cm⁻¹ s⁻² K⁻⁴, aT⁴/rho = g cm⁻¹ s⁻² K⁻⁴K⁴cm³g⁻¹ = cm² s⁻²
AU = 1.496e+13
mp = 1.6605402e-24 #g


#SANA/IMF PROB DISTRs.
alpha_p = -0.55
stdev_p=0.22
sigma_p=0
alpha_q = -0.10
stdev_q=0.58
sigma_q=0
alpha_m = -2.35


#useful symbols 
# ⭐ 🌀 ⚡ 💧 🔥 ✨ ⚽ ⚾ 🥌 🔮 🪄 ♠ ♥ ♦ ♣ 🔔 🎵 🎤 📀 🗑 💣 🧲 ⚠ ⛔ 🚫 ☢ ☣
# 🔆 ⭕ ✔ ✅ ❌ 🅰 🫟 🅱 🅾 🔴 🟠 🟡 🟢 ⚪ 🟣 🔵 🟤 ⚫ 🏁
Lsun = L☼ 
Msun = M☼
Rsun = R☼

#if working with Tuples instead of Dict, this is the similar way to get values
#value = getproperty(SN, Symbol(entry))   # → 1.0

function round_to_nearest_multiple(value, multiple, round_to = RoundNearest)
        return Base._round_invstep(value, 1/multiple, round_to)
end

function tick_logger(min, max)

        logmin, logmax = floor(log10(min)), ceil(log10(max))
        logs = collect(range(logmin,logmax, step=1))
        ints = collect(range(1,10,step=1))
        return [i*10^l for l in logs for i in ints ]
end  

function maxxer(arr)
        return [maximum(arr[1:i]) for i in 1:length(arr)]
end

# Function for interpolation with specified out-of-bounds handling
function interpolator(x, y; method=:linear, extrapolation=:flat, extrap_const=0)
        # Select interpolation method (only :linear and :constant for gridded interpolation)
        itp_method = method == :linear   ? Gridded(Linear()) :
                        method == :constant ? Gridded(Constant()) :
                        throw(ArgumentError("Invalid interpolation method for non-uniform grid"))

        # Create gridded interpolation object
        itp = interpolate((x,), y, itp_method)

        # Select extrapolation method
        if extrapolation == :flat
                extrap = extrapolate(itp, Flat())  # Constant value for out-of-bounds
        elseif extrapolation == :periodic
                extrap = extrapolate(itp, Periodic())      # Periodic extrapolation
        #elseif extrapolation == :line
        #        extrap = extrapolate(itp, Line())          # Linear extrapolation
        elseif extrapolation == :throw
                extrap = extrapolate(itp, NaN)          # Linear extrapolation
        elseif extrapolation == :constant
                extrap = extrapolate(itp, extrap_const)          # Linear extrapolation
        else
                throw(ArgumentError("Invalid extrapolation method"))
        end

    # Return the interpolation and extrapolation function
    return extrap
end


def_format_width = "%40"
function assigner(value; fmt=def_format_width, typ=".6e", quality_check = false)
        if isnan(value)
                return Dict("value" => NaN, "format" => fmt*typ)
        end 

        if quality_check
                if abs(value) >= 1e16 && 'f' in typ
                value = NaN 
                typ = ".6e"
                end 
        end 
        return Dict("value" => value, "format" => fmt*typ)
end 


function find_nearest(array, value; seek_ini=1, seek_end=-1, tol = 0.5, out_tol_value = 0, tol_scale = "abs")
        n = abs.(array[seek_ini : (seek_end == -1 ? end : seek_end)] .- value)
        idx = argmin(n)
        if tol_scale == "rel" && abs(n[idx]/value-1) >= tol
            #@printf("%.3e  vs %.3e -- %.3e vs %.3e\n", n[idx], value, n[idx]/value-1, tol)
            return out_tol_value
        end 
        if tol_scale == "abs" && (n[idx]-value) >= tol
            #@printf("found %.3e  value %.3e -- diff %.3e tol %.3e --- out %.3e\n", n[idx], value, n[idx]-value, tol, out_tol_value)
            return out_tol_value
        end 
    
        return (seek_ini-1) + idx
    end
    
function find_nearest_true(array, value; seek_ini=1, seek_end=-1, tol = 0.5, out_tol_value = 0, tol_scale = "abs")
        n = abs.(array[seek_ini : (seek_end == -1 ? end : seek_end)] .- value)
        idx = argmin(n)
        if tol_scale == "rel" && abs(n[idx]/value-1) >= tol
                #@printf("%.3e  vs %.3e -- %.3e vs %.3e\n", n[idx], value, n[idx]/value-1, tol)
                return out_tol_value
        end 
        if tol_scale == "abs" && (n[idx]-value) >= tol
                #@printf("found %.3e  value %.3e -- diff %.3e tol %.3e --- out %.3e\n", n[idx], value, n[idx]-value, tol, out_tol_value)
                return out_tol_value
        end 

        return seek_ini + idx
end
   
function kepler_law(; m, P, a)
 #       isnothing(m) 

 #       isnothing(P) 

        if isnothing(a) 
                return (cbrt.((P .* day).^ 2 ./ (4 .* pi .^ 2) .* G .* m .* Msun ) ./ Rsun)
        elseif isnothing(P) 
                return  sqrt( (a*Rsun)^3 * 4*pi^2 / (G*m*Msun))/day
        end 

        throw(ErrorException)
end

function eval_RL(m::AbstractFloat, m_companion::AbstractFloat, which::String, val::AbstractFloat) #m1 is the star for which you want to eval the RL
    #UNITS ARE IN MSUN, RSUN and DAY
    P = nothing 
    a = nothing 
    (which in ["P", "Period", "period"])       && (P = val) 
    (which in ["a", "separation", "distance"]) && (a = val)
    (P===nothing && a===nothing) && throw(ErrorException)
    if a===nothing
        a = cbrt.((P .* day).^ 2 ./ (4 .* pi .^ 2) .* G .* (m .+ m_companion) .* Msun ) ./ Rsun
    end
    qq = (m ./ m_companion) .^ (1/3)
    return  a .* 0.49 .* qq .* qq ./ (0.6 .* qq .* qq + log.(1 .+qq))
end

function eval_RLout(q)
        s = 49.4 ./ (12.2 .+ q .^ 0.208)
        return 1 .+ 2.74 ./ ( 1 .+ ((log.(q) .+ 1.02) ./ s) .^ 2 ) .* 1 ./ (7.13 .+ q .^ (-0.386))
end

function to_key(val)
        return @sprintf("%.2f", val)
    end
    


function read_table(filename)
        return CSV.read(filename, header=1, DataFrame,delim=' ',ignorerepeated=true)
end 


function evolutionary_checkpoints(model, ix_end; from_ZAMS=true)
        h1  = model.center_h1
        he4 = model.center_he4
        c12 = model.center_c12
        ZAMS=0
        TAMS=0
        end_Hburn = 0
        H_burn = Dict("99%"=>0, "0.50"=>0, "0.25"=>0)
        He_burn = Dict("0.75"=>0, "0.50"=>0, "0.25"=>0)
        ini_Heburn = 0 
        end_Heburn = 0
        ini_Cburn = 0
        end_Cburn = 0

        he4_endHburn = -1
        c12_endHeburn = -1

        end_burn_thresh = 1e-2 
        end_burn_thresh_2 = 1e-4
        ini_burn_thersh = 0.99

        for j in range(1, ix_end)
                if from_ZAMS
                        if (H_burn["99%"]==0) && h1[j] <= h1[1]*0.99
                                H_burn["99%"] = j 
                        end 
                        if (ZAMS == 0) && h1[j] <= 0.701
                                ZAMS = j 
                        end 
                        (ZAMS==0) && continue
                        if H_burn["0.50"]==0 && h1[j] <= 0.50
                                H_burn["0.50"]=j 
                        end
                        (H_burn["0.50"]==0) && continue
                        if H_burn["0.25"]==0 && h1[j] <= 0.25
                                H_burn["0.25"]=j 
                        end
                        (H_burn["0.25"]==0) && continue
                        if (TAMS==0) && h1[j] < end_burn_thresh
                                TAMS = j 
                        end 
                        if (end_Hburn == 0) && h1[j] < end_burn_thresh
                                end_Hburn = j
                        end 
                        (end_Hburn==0) && continue
                        if (he4_endHburn == -1) && h1[j] < end_burn_thresh_2
                                he4_endHburn = he4[j]
                        end 
                        (he4_endHburn>0) ? nothing : continue
                else
                        he4_endHburn =  he4[1]
                end 

                if (ini_Heburn == 0) && he4[j] < ini_burn_thersh*he4_endHburn
                        ini_Heburn = j
                end 
                (ini_Heburn==0) && continue
                if He_burn["0.75"]==0 && he4[j] <= 0.75
                        He_burn["0.75"]=j 
                end
                (He_burn["0.75"]==0) && continue
                if He_burn["0.50"]==0 && he4[j] <= 0.50
                        He_burn["0.50"]=j 
                end
                (He_burn["0.50"]==0) && continue
                if He_burn["0.25"]==0 && he4[j] <= 0.25
                        He_burn["0.25"]=j 
                end
                (He_burn["0.25"]==0) && continue
                if (end_Heburn == 0) && he4[j] < end_burn_thresh_2
                        end_Heburn = j
                end 
                end_Heburn>0 ? nothing : continue
                if (c12_endHeburn == -1) && he4[j] < end_burn_thresh_2
                        c12_endHeburn = c12[j]
                end 
                c12_endHeburn>0 ? nothing : continue
                if (ini_Cburn == 0) && c12[j] < ini_burn_thersh*c12_endHeburn
                        ini_Cburn = j
                end 
                ini_Cburn>0 ? nothing : continue
                if (end_Cburn == 0) && c12[j] < end_burn_thresh
                        end_Cburn = j
                end 
                end_Cburn>0 ? nothing : continue
        end 
        return ZAMS, H_burn, TAMS, end_Hburn, ini_Heburn, He_burn, end_Heburn, ini_Cburn, end_Cburn
end 

function str_significant(x::Float64, sigdig::Int)
        (x == 0) && (return "0")
        x = round(x, sigdigits=sigdig)
        n = length(@sprintf("%d", abs(x)))              # length of the integer part
        if (x ≤ -1 || x ≥ 1)
            decimals = max(sigdig - n, 0)               # 'sig - n' decimals needed 
        else
            Nzeros = ceil(Int, -log10(abs(x))) - 1      # No. zeros after decimal point before first number
            decimals = sigdig + Nzeros
        end 
        return @sprintf("%.*f", decimals, x)
    end
    

function fill_nan_linear_x(A, x_unique)
        A_filled = copy(A)
        for j in 1:size(A, 2)
                for i in 1:size(A, 1)
                if isnan(A[i, j])
                        # Search for neighbors in x-direction (i-axis)
                        left = findlast(!isnan, A[1:i-1, j])
                        right = findfirst(!isnan, A[i+1:end, j])
                        right = isnothing(right) ? nothing : right + i  # adjust index if found

                        if !isnothing(left) && !isnothing(right)
                        # Linear interpolation
                        x1, x2 = x_unique[left], x_unique[right]
                        z1, z2 = A[left, j], A[right, j]
                        x_missing = x_unique[i]
                        A_filled[i, j] = z1 + (z2 - z1) * (x_missing - x1) / (x2 - x1)
                        elseif !isnothing(left)
                        A_filled[i, j] = A[left, j]  # Nearest on the left
                        elseif !isnothing(right)
                        A_filled[i, j] = A[right, j]  # Nearest on the right
                        end
                end
                end
        end
        return A_filled
end














"""
    get_nearest_category(x, y, xq, yq)

Given arrays `x`, `y`, and  `z`, returns the index
at the nearest `(x, y)` point to the query point `(xq, yq)`.

Arguments:
- `x`, `y` : Arrays of coordinates of the known data points
- `xq`, `yq` : Query point

Returns:
- `index` of the nearest data point
"""
function get_nearest_category_weighted_ix(x, y,   xq, yq; xweight=1.0, yweight=1.0)
        @assert length(x) == length(y) 
    
        # Scale data according to weights
        xs = x .* xweight
        ys = y .* yweight
    
        data = hcat(xs, ys)'
        tree = KDTree(data)
    
        # Scale the query point too
        xq_scaled = xq * xweight
        yq_scaled = yq * yweight
    
        idxs, _ = knn(tree, [xq_scaled, yq_scaled], 1)
        return idxs[1]
end
    
killnan(x::Int; replacewith=0) =  isnan.(x) ? replacewith : x
killnan(x::Float64; replacewith=-Inf) =  isnan.(x) ? replacewith : x
killnan(x::Vector{Float64}; replacewith=-Inf) =  [isnan(x_n) ? replacewith : x_n for x_n in x ]
nanmax(a1,a2) =max.(killnan(a1), killnan(a2))
 

struct explosion_properties 
        M_ni::Float32
        E_exp::Float64 
        M_remnant_g::Float16 
        M_remnant_b::Float16
        v_kick::Float16
        function explosion_properties(M_ni::AbstractFloat, E_exp::Float64, M_remnant_g::AbstractFloat, M_remnant_b::AbstractFloat, v_kick::AbstractFloat)
                if any(isnan, (M_ni, E_exp, M_remnant_g, M_remnant_b, v_kick))
                    error("All values must be positive. Got: M_ni = $M_ni, E_exp = $E_exp, M_remnant_g = $M_remnant_g, M_remnant_b = $M_remnant_b, v_kick=$v_kick")
                end
                new(M_ni, E_exp, M_remnant_g, M_remnant_b, v_kick)
        end
        function explosion_properties(M)
                new(0., 0., M, M, 0)
        end
end

struct output_run_data
        M_end::Float16
        M_he::Float16
        M_co::Float16
        Xc::Float16
        deltaM_C::Float16
        M_ni::Float32
        E_exp::Float64 
        M_remnant_g::Float16 
        M_remnant_b::Float16
        v_kick::Float16
        solver_string::String
        function output_run_data(M_end, M_he, M_co, Xc, deltaM_C, exp_output; solver_string="X")
                if any(isnan, (M_end, M_he, M_co, Xc, deltaM_C))
                    error("All values must be positive. Got: M_end = $M_end, M_he = $M_he, M_co = $M_co, Xc = $Xc, deltaM_C = $deltaM_C")
                end
                
                new(M_end, M_he, M_co, Xc, deltaM_C, exp_output.M_ni, exp_output.E_exp,exp_output.M_remnant_g, exp_output.M_remnant_b, exp_output.v_kick, solver_string)
        end
        function output_run_data(M_end; solver_string="X")
                if any(isnan, (M_end))
                    error("All values must be positive. Got: M_end = $M_end")
                end
                new(M_end, M_end, M_end, 0, 0, 0, 0, M_end, M_end, 0, solver_string)
        end
        function output_run_data( )
                new(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, "X")
        end

end




function remove_nan(value_to_check_for_nan, value_to_replace_the_nan_with)

        isnan(value_to_check_for_nan) ? (return value_to_replace_the_nan_with) : (return value_to_check_for_nan)
    
end






print("done!\n")

    