function does_it_CC(;Mhe::Float64=-1., Mco::Float64=-1., Mconvhe::Float64=-1.,
                     Mhe_crit::Float64=2.49, 
                     Mco_crit::Float64=1.43, #1.435,
                     Mconvhe_crit::Float64=1.30)
    Mhe_criterion = false 
    Mco_criterion = false 
    Mconvhe_criterion = false 
    (Mhe     > 0) && (Mhe_criterion     = Mhe     > Mhe_crit    )
    (Mco     > 0) && (Mco_criterion     = Mco     > Mco_crit    )
    (Mconvhe > 0) && (Mconvhe_criterion = Mconvhe > Mconvhe_crit)
    
    return (Mhe_criterion || Mco_criterion || Mconvhe_criterion)
    
end 
print("\nBuilding Explosion criteria...")
print("PS20...")

PS2020_explode =  CSV.read(PS20_dir*"PS2020_explosion.dat", header=1, DataFrame,delim=' ',ignorerepeated=true)
PS2020_explode[!, "value"] .= "explode"
PS2020_implode =  CSV.read(PS20_dir*"PS2020_implosion.dat", header=1, DataFrame,delim=' ',ignorerepeated=true)
PS2020_implode[!, "value"] .= "implode"
PS2020 = vcat(PS2020_explode, PS2020_implode)
PS2020_extremes = Dict("Mco" => [minimum(PS2020.Mco), maximum(PS2020.Mco)], "Xc" => [minimum(PS2020.Xc), maximum(PS2020.Xc)]) 
function explosion_interpolation_metric(Mco, Xc; Mco_ref=PS2020.Mco, Xc_ref=PS2020.Xc)

    distance_x = (Mco .- Mco_ref) ./ (PS2020_extremes["Mco"][2]-PS2020_extremes["Mco"][1])
    distance_y = ( Xc .-  Xc_ref) ./ ( PS2020_extremes["Xc"][2]- PS2020_extremes["Xc"][1])

    return sqrt.( distance_x .^ 2 .+ distance_y .^ 2)
end 

function find_nearest_2D_explosion(Mco, Xc, Mco_ref=PS2020.Mco, Xc_ref=PS2020.Xc)
    Mco >= 10 && return "implode"
    distance_map = explosion_interpolation_metric(Mco, Xc; Mco_ref=PS2020.Mco, Xc_ref=PS2020.Xc)
    ix = argmin(distance_map)
    #@printf("point (%.3f,%.3f) mapped to (%.3f,%.3f)", Mco, Xc, Mco_ref[ix], Xc_ref[ix])
    #@printf(", which is found to %s", PS2020.value[ix] )
    return PS2020.value[ix] 
end


print("detailed (M16/Ertl)...")

fexp_single(M16_param) = CSV.read(CC_dir * "single_star/EXP_PROP_"* M16_param*".data",  
            DataFrame, delim=' ',ignorerepeated=true, silencewarnings=true)
fexp_HeS(M16_param)    = CSV.read(CC_dir*"HeS/EXP_PROP_"* M16_param*".data",  
          DataFrame, delim=' ',ignorerepeated=true, silencewarnings=true)
fexp_DRAD    = CSV.read(data_dir*"DRAD/EXP_PROP_AD23.data",  
          DataFrame, delim=' ',ignorerepeated=true, silencewarnings=true)

function corr_file_reducer_CC(corr_file, exp_file; mass_column = "logM")
    corr_file_exp = copy(corr_file)

    for ix in range(1,length(corr_file[!, mass_column]))
        logM = corr_file[!, mass_column][ix] 
        if ! (logM in exp_file.logM)
            #@printf("logM = %.3f -> is in exp? %s", logM, logM in exp_file.logM)
            #@printf(" -> removing from CORR")
            (mass_column == "logM") && (corr_file_exp = filter(row -> (row.logM in exp_file.logM), corr_file ))
            (mass_column == "M")    && (corr_file_exp = filter(row -> (row.M in exp_file.logM), corr_file ))

            #@printf("\n")
        end
    end
    # if (mass_column == "logM") 
    #     @printf("lenght corr %4d - exp %4d\n", length(corr_file_exp.logM), length(exp_file.logM))
    #     print("... ", corr_file_exp.logM == exp_file.logM, "\n")
    #     for logm in corr_file_exp.logM
    #         print(logm, ", ")
    #     end
    #     print("\n")
    #     for logm in exp_file.logM
    #         print(logm, ", ")
    #     end
    # end
    # (mass_column == "M")    &&  print("... ", corr_file_exp.M    == exp_file.logM, "\n")

    corr_file_exp[!, "E_exp"] = exp_file.E_exp  
    corr_file_exp[!, "M_ej"] = exp_file.M_ej  
    corr_file_exp[!, "xi25"] = exp_file.xi25  
    corr_file_exp[!, "mu4"] = exp_file.mu4  
    corr_file_exp[!, "M4"] = exp_file.M4  
    corr_file_exp[!, "Mni"] = exp_file[!, "M(Ni56)"]  
    corr_file_exp[!, "M_rem_b"] = exp_file.M_star .- exp_file.M_ej
    corr_file_exp[!, "M_rem_g"] = exp_file.M_grav
    corr_file_exp[!, "v_kick"] = exp_file.v_kick

    return corr_file_exp
end 
exp_corr_S = Dict() 
exp_corr_He = Dict()
M16_parameters = ["M16", "S21", "AD23"]
for parameters in M16_parameters 
    exp_corr_S[parameters]  = corr_file_reducer_CC(f_singleCC,        fexp_single(parameters); mass_column = "logM")
    exp_corr_He[parameters] = corr_file_reducer_CC(f_HeSCC,           fexp_HeS(parameters);    mass_column = "logM")
end 




print("...M imported CC grids, explodability files and interpolants\n")
wie_warning = 0
# Mhe, Mco, Xc, Mmax, when, criterion, reference_grid; debug=false
function will_it_explode(;Mhe::Float64=0., Mco::Float64=0., Xc::Float64=0., Mmax::Float64=0., when::String="", 
                            criterion::String="", reference_grid::String="", 
                            debug = false, Mconv_HeBurn = NaN)
    #this function seeks the model in the reference grid with the closest He and CO core masses. 
    #if the chosen grid is the MW-ZAMS grid, then the maximum initial mass must be also included
    #so as to break the degeneracy between low and high-mass models that end up with the same mass
    
    if (Mhe<0 || Mco <0 || Xc <0 || Mmax<0 || length(when)<=0 || length(criterion)<=0 || length(reference_grid)<=0)
        @printf("Input values\nMhe=%5.3f\nMco=%5.3f\nXc=%5.3f\nMmax=%5.3f\nwhen=%s\ncriterion=%s\nreference_grid=%s",
                    Mhe, Mco, Xc, Mmax, when, criterion, reference_grid)
        throw(ErrorException("Some variable you provided is not ok (will_it_explode). Check them again. Terminating!"))

    end

    global wie_warning
    
    if ! (occursin("MM", criterion) || occursin("Ertl", criterion) || criterion in ["PS20", "X"] || occursin("xi", criterion) || occursin("comp", criterion))
        throw(ErrorException("EXPLOSION CRITERION ILL-DEFINED. CHECK WHAT YOU ENTERED\n"))
    end 

    M16_param = occursin("MM", criterion) ? criterion[4:end] : "M16"  
    fexp = nothing 
    extrapolate_lowM = false
    trigger_extrapolate_lowM= false
    explosion_params = nothing 
    ix = 0
    cut_range = nothing
    # print("read criterion ", criterion, "\n")
    if reference_grid == "MW-ZAMS"
        fexp = exp_corr_S[M16_param]
        logM = fexp.logM
    
        cut1 = find_nearest(fexp.logM, 1.59)
        cut2 = find_nearest(fexp.logM, 1.66)
        cut3 = find_nearest(fexp.logM, 1.80)
        range1 =    1:cut1
        range2 = cut1:cut2
        range3 = cut2:cut3
        range4 = cut3:length(fexp.logM)
        ranges = [range1, range2, range3, range4]
        if log10(Mmax) <= 1.59
            cut_range = range1
            ix_range = 1
        elseif 1.59 < log10(Mmax) <= 1.66
            cut_range = range2
            ix_range = 2
        elseif 1.66 < log10(Mmax) <= 1.80
            cut_range = range3
            ix_range= 3
        elseif  1.80 < log10(Mmax) 
            cut_range = range4
            ix_range = 4
        end
        debug && print(cut_range, " ", cut )
    elseif reference_grid == "MW-Stripped"
        fexp = exp_corr_He[M16_param]
        cut_range = 1:1:length(fexp.E_exp)
        ranges = [cut_range]
        ix_range=0
        extrapolate_lowM = true
    elseif reference_grid == "AguileraDena"
        fexp = exp_corr_AD
        cut_range = 1:1:length(fexp.E_exp)
        ranges = [cut_range]
        ix_range= 0
    end

    extra_threshold = 0.05
    try
        cut_i = first(cut_range)
        cut_f =  last(cut_range)

        if !isnan(Mco) && Mco >= maximum(fexp[!, "M_co_core_"*when])
            ix = find_nearest(fexp[!, "M_co_core_"*when], maximum(fexp[!, "M_co_core_"*when]))
            which = "M_co_core"
            x = Mco
        elseif !isnan(Mhe) && Mhe >= maximum(fexp[!, "M_he_core_"*when])
            ix = find_nearest(fexp[!, "M_he_core_"*when], maximum(fexp[!, "M_he_core_"*when]))
            which = "M_he_core_"
            x = Mhe
        elseif !isnan(Mco) && minimum(fexp[!, "M_co_core_"*when][cut_range]) <= Mco <= maximum(fexp[!, "M_co_core_"*when][cut_range])
            ix = find_nearest(fexp[!, "M_co_core_"*when], Mco, seek_ini=minimum(cut_range), seek_end=maximum(cut_range))
            which = "M_co_core"
            x = Mco
        elseif !isnan(Mhe) && minimum(fexp[!, "M_he_core_"*when][cut_range]) <=  Mhe <= maximum(fexp[!, "M_he_core_"*when][cut_range])
            ix = find_nearest(fexp[!, "M_he_core_"*when], Mhe, seek_ini=minimum(cut_range), seek_end=maximum(cut_range))
            which = "M_he_core"
            x = Mhe
        elseif !isnan(Mco) && Mco < minimum(fexp[!, "M_co_core_"*when][cut_range] )
            which = "M_co_core"
            x = Mco
            if cut_i != 1 
                # print("UNDERMASSIVE M_CO FOR THE RANGE")
                # throw(ErrorException)
                ix = find_nearest(fexp[!, "M_co_core_"*when], Mco)
                (ix == 0) && (ix = 1)
                # @printf("WARNING! Explosion association was unpredictable.\n")
                wie_warning += 1
            else 
                ix = cut_i 
            end
        elseif !isnan(Mhe) && Mhe < minimum(fexp[!, "M_he_core_"*when][cut_range])
            which = "M_he_core"
            x = Mhe
            if cut_i != 1 
                print("UNDERMASSIVE M_HE FOR THE RANGE")
                throw(ErrorException)
            else 
                ix = cut_i 
            end
        elseif !isnan(Mco) && Mco > maximum(fexp[!, "M_co_core_"*when][cut_range] )
            which = "M_co_core"
            x = Mco
            if cut_f != length(fexp.E_exp)
                # print("OVERMASSIVE M_CO FOR THE RANGE")
                # throw(ErrorException)
                ix = find_nearest(fexp[!, "M_co_core_"*when], Mco)
                # @printf("WARNING! Explosion association was unpredictable.\n")
                (ix == 0) && (ix = 1)
                wie_warning += 1
            else 
                ix = cut_f
            end
        elseif !isnan(Mhe) && Mhe > maximum(fexp[!, "M_he_core_"*when][cut_range])
            which = "M_he_core"
            x = Mhe
            if cut_f != length(fexp.E_exp)
                print("OVERMASSIVE M_HE FOR THE RANGE")
                throw(ErrorException)
            else 
                ix = cut_f
            end
        else
            #undermassive!
            throw(ErrorException("Keep Working"))
            trigger_extrapolate_lowM=true
            ix = 1
        end 
        ix == 0 && throw(ErrorException)
    catch e
        @printf("%s\n!!! ERROR in finding the closest explosion model", e)
        @printf("ix = %d, Mco = %.2f, Mhe = %.2f, which = %s, when = %s,\n reference_grid = %s, extrapolate_lowM = %s\n",
        ix, Mco, Mhe, which, when, reference_grid, extrapolate_lowM)
        @printf("log10(Mmax)=%.2f - cut range : %.2f - %.2f\n",log10(Mmax), fexp.logM[first(cut_range)], fexp.logM[last(cut_range)] )
        @printf("limits: Mco/Mhe between = %.2f/%.2f %.2f/%.2f\n", 
                    fexp[!, "M_co_core_"*when][first(cut_range)], fexp[!, "M_he_core_"*when][first(cut_range)],
                    fexp[!, "M_co_core_"*when][last(cut_range)], fexp[!, "M_he_core_"*when][last(cut_range)],
                    )



        throw(ErrorException("Terminating"))
    end 

    E_exp = NaN
    if trigger_extrapolate_lowM

        E_exp = NaN
        if extrapolate_lowM
            if !isnan(Mco) && Mco>1.44
                E_exp = 1e50 + (fexp.E_exp[1]-1e50)/(fexp[!, "M_co_core_"*when][1]-1.44) * (Mco-1.44)
            elseif !isnan(Mhe) && Mhe > 2.49
                E_exp = 1e50 + (fexp.E_exp[1]-1e50)/(fexp[!, "M_he_core_"*when][1]-2.49) * (Mhe-2.49)
            elseif  !isnan(Mconv_HeBurn)
                E_exp = 1e50 
            else 
                throw(ErrorException)
            end
        else 
            E_exp = fexp.E_exp[1]
        end 

        explosion_params = explosion_properties(fexp.Mni[1], E_exp, fexp.M_rem_g[1], fexp.M_rem_b[1], fexp.v_kick[1])
    else 
        E_exp = fexp.E_exp[ix]
        Mni = fexp.Mni[ix]
        M_rem_g = fexp.M_rem_g[ix]
        M_rem_b = fexp.M_rem_b[ix]
        v_kick  = fexp.v_kick[ix]
        if E_exp <= 0 
            Mni = 1e-2
            M_rem_g = 1.4 
            M_rem_b = 1.4 
            v_kick = 1e99
        end
        explosion_params = explosion_properties(Mni, E_exp, M_rem_g, M_rem_b, v_kick, "from will_it_explode")

    end

    # @printf("reference_grid = %s, which %20s, when %20s, mrange %10s, x %5.3f\n",reference_grid,  which, when, Mrange, x)

    if occursin("MM", criterion) 
        debug && @printf("Eexp = %6.3e, Mco = %4.1f  (found[%3d]  %4.1f, diff = %5.2f)", fexp.E_exp[ix], Mco, ix, fexp[!, "M_co_core_"*when][ix], Mco-fexp[!, "M_co_core_"*when][ix])
        return (E_exp > 0, explosion_params) 
    elseif occursin("Ertl", criterion) 
        Ertl_calibration = criterion[6:end]
        ks = Dict(
        "S19.8_old" =>[0.274, 0.0470],
        "S19.8_new" =>[0.294, 0.0468],
        "W15_old"   =>[0.225, 0.0495], 
        "W15_new"   =>[0.225, 0.0521], 
        "W18_old"   =>[0.283, 0.0430], 
        "W18_new"   =>[0.283, 0.0438], 
        "W20_old"   =>[0.284, 0.0393],
        "W20_new"   =>[0.273, 0.0413],
        "N20_old"   =>[0.194, 0.0580],
        "N20_new"   =>[0.182, 0.0608]
        )
        k=ks[Ertl_calibration]
        mu4 = fexp.mu4[ix]
        M4 = fexp.M4[ix]
        x = mu4 * M4
        y = mu4
        ycrit = k[1] .* x .+ k[2]
        debug && @printf("Mco = %4.1f  (found[%3d]  %4.1f, diff = %5.2f) M4mu4=%7.3f  mu4 = %7.3f mu4_crit=%7.3f ", Mco, ix, fexp[!, "M_co_core_"*when][ix], Mco-fexp[!, "M_co_core_"*when][ix], x, y, ycrit)

        return (y .< ycrit, explosion_params )
    elseif criterion == "X"
        return (true, explosion_params) 
    elseif criterion == "PS20"
        if (isnan(Mco) || isnan(Xc))
            @printf(" recieved Mco = %5.3f & Xc = %5.3f\n", Mco, Xc)
            throw(ErrorException("FOR PS20, Mco and Xc SHOULDNT BE NANs"))
        else
            does_it_explode = ( find_nearest_2D_explosion(Mco, Xc) == "explode" )
            return (does_it_explode, explosion_params) 

        end
    elseif occursin("comp", criterion)

        threshold = parse(Float64, criterion[5:end])
        xi = fexp.xi25[ix]
        return (xi < threshold, explosion_params)

    end

    print("EXPLOSION CRITERION ILL-DEFINED. CHECK WHAT YOU ENTERED\n")
    throw(ErrorException)

end

