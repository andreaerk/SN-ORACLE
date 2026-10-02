

print("Reading & cleaning single-star CORR files...")
f_single = CSV.read(SG_dir*"CORR_GRID_single_150.dat", 
                     header=1, DataFrame,delim=' ',ignorerepeated=true)
deleteat!(f_single, findall( !>(0) , f_single.M_he_core_end_BURN_He) )

f_singleCC = CSV.read(CC_dir*"single_star/CORR_GRID_single_toCC.dat",  DataFrame, delim=' ',ignorerepeated=true)
f_HeSCC = CSV.read(CC_dir*"HeS/CORR_GRID_single_toCC.dat",  DataFrame, delim=' ',ignorerepeated=true)
print("done!\n")

function merge_corr_files(dataframe1, dataframe2)
    set = names(dataframe2)

    dataframe_output = copy(dataframe1[!,set  ])
    #global dataframe_tmp = copy(dataframe1[!,set  ])
    dataframe_tmp = copy(dataframe1[!,set  ])
    for ix in range(1,length(dataframe_output.logM))
        #print(f_single3)
        logM = dataframe_output.logM[ix] 
        if (logM in dataframe2.logM)
            #@printf("logM = %.3f -> is in exp? %s", logM, logM in exp_file.logM)
            #@printf(" -> removing from CORR")
            filter!(row -> !(row.logM in dataframe2.logM), dataframe_tmp )

            #@printf("\n")
        end
    end
    #print(f_single3)
    dataframe_output = vcat(dataframe_tmp, dataframe2)

    return dataframe_output
end

f_single2 =  merge_corr_files(f_single, f_singleCC)
            

function produce_logRl_m(HeS_file)
    f_R_HeS = CSV.read(HeS_file,
                            header=1, DataFrame,delim=' ',ignorerepeated=true) 
    ts_HeS = ["max", "1yr", "1kyr", "2kyr", "5kyr", "10kyr", "20kyr"]
    Mmxconv_HeS = f_R_HeS.Mconvmax_Heburn
    log_RL_m = Dict()
    for t in ts_HeS 
        log_RL_m[t] = interpolator(Mmxconv_HeS, log10.(f_R_HeS[!, "R"*t]), extrapolation=:flat)
    end 

    return log_RL_m
end

log_RL_m = produce_logRl_m(CC_dir*"/He_star_summary_Ercolino2025.dat")

function predictor_caseBC(f, model; which="1", t_treshold = "1yr")
    Mconv_HeBurn = f[!, "Mconv_max_Heburn_"*which][model]
    if isnan(Mconv_HeBurn)
        return false
    end 

    if ! isnan( f[!, "age_"*which*"_end_BURN_He_"*which][model])
        RL = f[!, "Rl_"*which*"_end_BURN_He_"*which][model]
    elseif ! isnan(f[!, "age_"*which*"_ini_BURN_He_"*which][model])
        RL = f[!, "Rl_"*which*"_end_"*which][model]
    else 
        print("Warning: model terminated too soon to estimate CaseBC RLOF - setting it to False")
        return false
    end

    Rmax = 10. .^ log_RL_m[t_treshold](Mconv_HeBurn)

    if RL >  Rmax #will not fill RLOF
        return false
    elseif RL <= Rmax #will fill RLOF
        return true
    else #WARNING! NAN
        #@printf("model = %5d: s:%1s,  M = %7.2f, Rl = %7.2f vs Rmax = %7.2f\n", model, which,  Mconv_HeBurn,  RL, Rmax)
        #print("Warning: model terminated too soon to estimate CaseBC RLOF - setting it to False")
        #throw(ErrorException("NaN value encountered for CaseBC predictor"))
        return false
    end
end 

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



print("Building merger Mi+f_acc->Mco, BSG interpolators: Retrieving data ...")
f24=CSV.read(S24_dir*"table_models_zenodo.csv", DataFrame)
deleteat!(f24, findall( ismissing , f24[!, " t_cc/Myr"]) )
function build_f24_interpolators(early_or_late, f=f24)

    f24_filter = early_or_late == "early" ? "e" : "l"
    to_use =  (" Case B-"*f24_filter) .== f24[!, " Case"] 

    singles = ( " Single" .== f[!, " Case"])
    Mco  =  convert.(Float64,  f[!, " M_CO/Msun"] )
    Mhe  =  convert.(Float64,  f[!, " M_He/Msun"] )
    Mini =  convert.(Float64, f[!, "M_ini/Msun"]) 
    Mend =  convert.(Float64, f[!, " M_final/Msun"]) 
    Xc  = f[!, " X_C"] 

    Mend_single = interpolator(Mini[singles], Mend[singles], extrapolation=:flat)
    Mco_single  = interpolator(Mini[singles], Mco[singles],  extrapolation=:flat)
    Mhe_single  = interpolator(Mini[singles], Mhe[singles],  extrapolation=:flat)
    Xc_single   = interpolator(Mini[singles], Xc[singles],   extrapolation=:flat)
    
    
    f_acc = [k == " " ? 0. : parse(Float64, k)  for k in f[!, " f_acc"]] 
    logTeff =  convert.(Float64, f[!, " log Teff_cc/K"]) 
    Mend_ratio = Mend ./Mend_single(Mini)
    Mco_ratio  = Mco  ./Mco_single(Mini)
    Mhe_ratio  = Mhe  ./Mhe_single(Mini)
    Xc_ratio   = Xc   ./Xc_single(Mini)

    x = Mini[to_use]
    y = f_acc[to_use]
    x_unique = sort(unique(x))
    y_unique = sort(unique(y))
    nx = length(x_unique)
    ny = length(y_unique)

    z_mf_mer   = Mend_ratio[to_use]
    z_mco_mer  = Mco_ratio[to_use]
    z_xc_mer   = Xc_ratio[to_use]
    z_Teff_mer = logTeff[to_use] 
    z_Mhe_mer  = Mhe_ratio[to_use]


    A_mf   = fill(NaN, nx, ny)  # fill missing values with NaN for now
    A_mco  = fill(NaN, nx, ny)  # fill missing values with NaN for now
    A_xc   = fill(NaN, nx, ny)  # fill missing values with NaN for now
    A_Teff = fill(NaN, nx, ny)  # fill missing values with NaN for now
    A_mhe  = fill(NaN, nx, ny)  # fill missing values with NaN for now

    x_idx = Dict(val => i for (i, val) in enumerate(x_unique))
    y_idx = Dict(val => i for (i, val) in enumerate(y_unique))
    for (xi, yi, zi) in zip(x, y, z_mco_mer)
        A_mco[x_idx[xi], y_idx[yi]] = zi
    end
    for (xi, yi, zi) in zip(x, y, z_xc_mer)
        A_xc[x_idx[xi], y_idx[yi]] = zi
    end
    for (xi, yi, zi) in zip(x, y, z_Teff_mer)
        A_Teff[x_idx[xi], y_idx[yi]] = zi
    end
    for (xi, yi, zi) in zip(x, y, z_Mhe_mer)
        A_mhe[x_idx[xi], y_idx[yi]] = zi
    end
    for (xi, yi, zi) in zip(x, y, z_mf_mer)
        A_mf[x_idx[xi], y_idx[yi]] = zi
    end
    
    A_filled_mco  = fill_nan_linear_x(A_mco, x_unique)
    A_filled_xc   = fill_nan_linear_x(A_xc, x_unique)
    A_filled_Teff = fill_nan_linear_x(A_Teff, x_unique)
    A_filled_mhe  = fill_nan_linear_x(A_mhe, x_unique)
    A_filled_mf   = fill_nan_linear_x(A_mf, x_unique)
    
    # Step 4: Create interpolator
    nodes = (x_unique, y_unique)
    mf_to_mco  = extrapolate(interpolate(nodes, A_filled_mco,  Gridded(Linear())), Flat())
    mf_to_xc   = extrapolate(interpolate(nodes, A_filled_xc,   Gridded(Linear())), Flat())
    mf_to_teff = extrapolate(interpolate(nodes, A_filled_Teff, Gridded(Linear())), Flat())
    mf_to_mhe  = extrapolate(interpolate(nodes, A_filled_mhe, Gridded(Linear())), Flat())
    mf_to_mf   = extrapolate(interpolate(nodes, A_filled_mf, Gridded(Linear())), Flat())
    

    return mf_to_mco, mf_to_xc, mf_to_teff, mf_to_mhe, mf_to_mf

end
print("building interpolators...")
mf_to_mco_e, mf_to_xc_e, mf_to_teff_e, mf_to_mhe_e, mf_to_mf_e = build_f24_interpolators("early")
mf_to_mco_l, mf_to_xc_l, mf_to_teff_l, mf_to_mhe_l, mf_to_mf_l = build_f24_interpolators("late")
print("done!\n")



print("Building merger/partial stripped outcome -> core He dep. Retrieving data ...")
fms = CSV.read(BG_dir*"Mcore_Menv_corr_new.dat",
                 DataFrame,delim=' ', ignorerepeated=true)
_logMhe=log10.(fms.Mhe)
_logMenv=log10.(fms.Menv)

counter_call_gefibH = 0

"""
Assume you have a star, be it one that crashed or merged, of which you know its info at core He ignition
Here, given a prebuilt database of partially-stripped model's evolution, find the model that fits closest,
given the same He-core mass and envelope mass, and spit out the outcome. 
"""
function get_end_from_ini_BURN_He(M_he_core, M_env, criterion; logMcore = _logMhe, logMenv=_logMenv, reference_grid="MW-ZAMS", caseC = false, solver_string = "", II_I_MH_threshold= 0.001, debug = false)
    
    global counter_call_gefibH += 1

    ix = get_nearest_category_weighted_ix(logMcore,logMenv, log10(M_he_core), log10(M_env); yweight=1/5)
    debug && @printf("index retrieved = %5d\n", ix)
    M_end = fms.Mend[ix]
    M_hecore_end = fms.Mhecore_end[ix]
    M_cocore_end = fms.Mcocore_end[ix]
    Xc =  fms.Xc[ix]
    Mmax = (reference_grid=="MW-ZAMS") ? M_he_core+M_env : M_he_core

    CC = does_it_CC(Mhe=M_hecore_end, Mco=M_cocore_end)
    debug && @printf("CC = %s, Mhecore = %5.3f, Mcocore = %5.3f\n", CC, M_hecore_end, M_cocore_end)
    (!CC) ? (return :WD, output_run_data(min(M_end,1.4); solver_string=solver_string) ) : nothing 

    explode, exp_par = will_it_explode(;   Mhe=M_hecore_end, Mco=M_cocore_end, 
                                                  Xc=Xc, Mmax=Mmax, 
                                                  when="end_BURN_He", criterion=criterion, 
                                                  reference_grid=reference_grid) 
   

    caseC && print("                  ^^^case C! IIn!\n")
    Mh = (M_end - M_hecore_end)
    Type=nothing
    if explode 
        ( caseC && Mh>=1) && (Type = :IIn)
        (!caseC && Mh>=1) && (Type = :IIP)
        (!caseC && II_I_MH_threshold<=Mh<1) && (Type = :IIb)
        (!caseC && Mh<II_I_MH_threshold  )   && (Type = :Ibc)
    else 
        Type = :BH
    end
    return Type, output_run_data(M_end, M_hecore_end, M_cocore_end, Xc, 0. , exp_par, solver_string=solver_string)

end


print("Building interpolation library...")
M_max_single = 10 .^ f_single.logM
M_max_singleCC = 10 .^ f_single2.logM
M_max_HeS = f_HeSCC.M_he_core_ini_BURN_He
_M_core = f_single.M_he_core_ini_BURN_He
M_core_diff = vcat([0],  diff(_M_core))
M_core_ix = [i for i in range(1,findall( !>(0) , M_core_diff)[2]-1)]
# deleteat!(M_core_ix, findall( !>(0) , M_core_diff))

f_age_MS   = interpolator(M_max_single, f_single.age_end_BURN_H,   extrapolation=:flat)
f_Rmax_MS  = interpolator(M_max_single, f_single.Rmax_end_BURN_H,    extrapolation=:flat)
f_age_pMS  = interpolator(M_max_single, f_single.age_ini_BURN_He,  extrapolation=:flat)
f_Rmax_pMS = interpolator(M_max_single, f_single.Rmax_ini_BURN_He,   extrapolation=:flat)
f_age_HeB = interpolator(M_max_single, f_single.age_end_BURN_He,  extrapolation=:flat)
f_Rmax_HeB= interpolator(M_max_single, f_single.Rmax_end_BURN_He,   extrapolation=:flat)
f_age_pHeB = interpolator(M_max_single, f_single.age_ini_BURN_C,  extrapolation=:flat)
f_Rmax_pHeB= interpolator(M_max_single, f_single.Rmax_ini_BURN_C,   extrapolation=:flat)
f_age_CB   = interpolator(M_max_single, f_single.age_end,   extrapolation=:flat)
f_Rmax_CB  = interpolator(M_max_single, f_single.Rmax_end,    extrapolation=:flat)
t_preCC_pMS_M = interpolator(M_max_single, f_single.age_end - f_single.age_ini_BURN_He, extrapolation=:flat)
t_preCC_pMS_Mcore  = interpolator(_M_core[M_core_ix], f_single.age_end[M_core_ix] - f_single.age_ini_BURN_He[M_core_ix], extrapolation=:flat)
# Mconvcore_pMS_Mcore  = interpolator(_M_core[M_core_ix], f_singleCC.Mconv_max_He[M_core_ix], extrapolation=:flat)
t_preCC_pHeB_Mcore = interpolator(_M_core[M_core_ix], f_single.age_end[M_core_ix] .- f_single.age_end_BURN_He[M_core_ix], extrapolation=:flat)


interpol_single_Mmax=Dict()
interpol_single_Mmax["M_end"]        = interpolator(M_max_singleCC, f_single2.M_end,   extrapolation=:flat)
interpol_single_Mmax["Mconv_max_He"] = interpolator(M_max_singleCC, f_single2.Mconv_max_He,   extrapolation=:flat)
interpol_single_Mmax["M_env_end"]    = interpolator(M_max_singleCC, f_single2.M_end-f_single2.M_he_core_end,   extrapolation=:flat)
interpol_single_Mmax["M_h_end"]      = interpolator(M_max_singleCC, f_single2.M_tot_h_end,   extrapolation=:flat)
interpol_single_Mmax["M_hecore_end"] = interpolator(M_max_singleCC, f_single2.M_he_core_end, extrapolation=:flat)
interpol_single_Mmax["M_cocore_end"] = interpolator(M_max_singleCC, nanmax.(f_single2.M_co_core_ini_BURN_C,f_single2.M_co_core_end), extrapolation=:flat)
interpol_single_Mmax["X_C"]          = interpolator(M_max_singleCC, f_single2.Cc_end_BURN_He, extrapolation=:flat)
interpol_single_Mmax["M_hecore_pMS"] = interpolator(M_max_singleCC, f_single2.M_he_core_ini_BURN_He, extrapolation=:flat)
interpol_single_Mmax["M_pMS"] = interpolator(M_max_singleCC, f_single2.M_end_BURN_H, extrapolation=:flat)
interpol_single_Mmax["logTeff"] = interpolator(M_max_singleCC, f_single2.logT_end, extrapolation=:flat)

interpol_HeS_Mmax=Dict()
interpol_HeS_Mmax["M_end"]        = interpolator(M_max_HeS, f_HeSCC.M_end,   extrapolation=:flat)
interpol_HeS_Mmax["Mconv_max_He"] = interpolator(M_max_HeS, f_HeSCC.Mconv_max_He,   extrapolation=:flat)
interpol_HeS_Mmax["M_hecore_end"] = interpolator(M_max_HeS, f_HeSCC.M_he_core_end, extrapolation=:flat)
interpol_HeS_Mmax["M_cocore_end"] = interpolator(M_max_HeS, f_HeSCC.M_co_core_end_BURN_C, extrapolation=:flat)
interpol_HeS_Mmax["X_C"]          = interpolator(M_max_HeS, f_HeSCC.Cc_end_BURN_He, extrapolation=:flat)

interpol_HeS_Mhedep=Dict()
ix1 = [i for i in findall( (!=(1.00))   , f_HeSCC.logM)]
ix2 = [i for i in findall( (!=(1.06))   , f_HeSCC.logM)]
ix = intersect(ix1, ix2)
interpol_HeS_Mhedep["M_cocore_end"] = interpolator(f_HeSCC.M_he_core_end_BURN_He[ix], f_HeSCC.M_co_core_end_BURN_C[ix], extrapolation=:flat)
interpol_HeS_Mhedep["X_C"]          = interpolator(f_HeSCC.M_he_core_end_BURN_He[ix], f_HeSCC.Cc_end_BURN_He[ix]/0.99, extrapolation=:flat)

interpol_single_M_pMS=Dict()
interpol_single_M_pMS["M_end"] = interpolator(f_single2.M_ini_BURN_He[M_core_ix], f_single2.M_end[M_core_ix],   extrapolation=:flat)
interpol_single_M_pMS["M_env_end"] = interpolator(f_single2.M_ini_BURN_He[M_core_ix], f_single2.M_end[M_core_ix]-f_single2.M_he_core_end[M_core_ix],   extrapolation=:flat)

interpol_single_Mcore_pMS=Dict()
interpol_single_Mcore_pMS["M_end"] = interpolator(_M_core[M_core_ix], f_single2.M_end[M_core_ix],   extrapolation=:flat)
interpol_single_Mcore_pMS["M_pMS"] = interpolator(_M_core[M_core_ix], f_single2.M_ini_BURN_He[M_core_ix],   extrapolation=:flat)
interpol_single_Mcore_pMS["M_Max"] = interpolator(_M_core[M_core_ix], M_max_singleCC[M_core_ix],   extrapolation=:flat)
interpol_single_Mcore_pMS["M_env_end"] = interpolator(_M_core[M_core_ix], f_single2.M_end[M_core_ix]-f_single2.M_he_core_end[M_core_ix],   extrapolation=:flat)



print("done!\nBuilding Explosion criteria...")

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
print("...PS20")



fexp_single(M16_param) = CSV.read(CC_dir * "single_star/EXP_PROP_"* M16_param*".data",  
            DataFrame, delim=' ',ignorerepeated=true)
fexp_HeS(M16_param)    = CSV.read(CC_dir*"HeS/EXP_PROP_"* M16_param*".data",  
          DataFrame, delim=' ',ignorerepeated=true)

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




print("...M16-E16 - imported CC grids, explodability files and interpolants\n")
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
    # retriever = nothing    
    # Mrange = nothing

    # logMmax = log10.(Mmax)
    
    # if reference_grid == "MW-ZAMS"
    #     fexp = exp_corr_S
    #     (         logMmax .<= 1.43) && (Mrange = "lowM")
    #     (2.00 .>= logMmax .>= 1.78) && (Mrange = "higM")
    #     (1.43 .<  logMmax .<  1.78) && (Mrange = "medM")
    #     (         logMmax .>  2.00) && (Mrange = "higM2")
    #     retriever = interpolator_singlestar_exp_params
    #     isnothing(Mrange) && (@printf("MHMHMH read %.3f", logMmax); throw(ErrorException))
    # elseif reference_grid == "MW-Stripped"
    #     fexp = exp_corr_He
    #      (        logMmax .>=  1.06)  && (Mrange = "higM")
    #      (0.98 .< logMmax .<   1.06)  && (Mrange = "medM")
    #      (        logMmax .<=  0.98)  && (Mrange = "lowM")
    #     retriever = interpolator_singlehestar_exp_params
    #     extrapolate_lowM = true 

    # else 
    #     throw(ErrorException("!GOT WRONG REFERENCE GRID!"))    
    # end

    # x = nothing 
    # which = nothing
    # if !isnan(Mco)
    #     which = "M_co_core"
    #     x = Mco
    # elseif !isnan(Mhe)
    #     which = "M_he_core"
    #     x = Mhe
    # else 
    #     #WHY ARE YOU RECEIVING NANS?
    #     throw(ErrorException("BOTH MCO and MHE ARE NANS! THAT SHOULD NOT HAPPEN!"))
    # end

    # trigger_extrapolate_lowM = (x <= fexp[!, which*"_"*when][1])

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
        explosion_params = explosion_properties(Mni, E_exp, M_rem_g, M_rem_b, v_kick)

        # E_exp     = retriever[which]["E_exp"][when][Mrange](x)
        # Mni       = retriever[which]["Mni"][when][Mrange](x)
        # M_rem_g   = retriever[which]["M_rem_g"][when][Mrange](x)
        # M_rem_b   = retriever[which]["M_rem_b"][when][Mrange](x)

        # try 
        #     explosion_params = explosion_properties(Mni, E_exp, M_rem_g, M_rem_b)
        # catch 
        #     print(reference_grid, which, when, Mrange, x, " ", logMmax, "\n")
        #     throw(ErrorException)
        # end 
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



function SINGLE_SN(M, criterion, kind; f_acc = 0, Mi=nothing, early_or_late = nothing, caseC = false, solver_string="")

    SN_explosion = true
    H = -1
    output = nothing 
    BSG = nothing 

    if kind == "stripped"
        reference_grid = "MW-Stripped"

        M_hecore      = interpol_HeS_Mmax["M_hecore_end"](M) 
        H             = 0.
        M_end         = interpol_HeS_Mmax["M_end"](M)
        M_COcore      = interpol_HeS_Mmax["M_cocore_end"](M)
        X_c_endburnHe = interpol_HeS_Mmax["X_C"](M)

        CC = does_it_CC(Mhe=M_hecore, Mco=M_COcore)

        (!CC) ? (return :WD,  output_run_data(min(M_end,1.4); solver_string=solver_string) ) : nothing  
        SN_explosion, exp_param = will_it_explode(;Mhe=M_hecore, Mco=M_COcore, 
                                                  Xc= X_c_endburnHe, Mmax=M, 
                                                  when="end_BURN_C", criterion=criterion, 
                                                  reference_grid=reference_grid)
        output = output_run_data(M_end, M_hecore, M_COcore, X_c_endburnHe, caseC ? Inf : 0, exp_param; solver_string=solver_string)
        BSG = false 
    elseif kind == "single"
        reference_grid = "MW-ZAMS"

        M_hecore      = interpol_single_Mmax["M_hecore_end"](M) 
        H             = interpol_single_Mmax["M_env_end"](M)
        M_end         = interpol_single_Mmax["M_end"](M)
        M_COcore      = interpol_single_Mmax["M_cocore_end"](M)
        X_c_endburnHe = interpol_single_Mmax["X_C"](M)

        CC = does_it_CC(Mhe=M_hecore, Mco=M_COcore)
        (!CC) ? (return :WD,  output_run_data(min(M_end,1.4); solver_string=solver_string) ) : nothing  
        SN_explosion, exp_param = will_it_explode(;Mhe=M_hecore, Mco=M_COcore, 
                                                  Xc= X_c_endburnHe, Mmax=M, 
                                                  when="end_BURN_C", criterion=criterion, 
                                                  reference_grid=reference_grid) 
        output = output_run_data(M_end, M_hecore, M_COcore, X_c_endburnHe, caseC ? Inf : -1, exp_param; solver_string=solver_string)
        BSG = interpol_single_Mmax["logTeff"](M) > 3.9
    elseif kind == "accretor" 
        reference_grid = "MW-ZAMS"

        mf_to_mco  = early_or_late == "early" ? mf_to_mco_e  : mf_to_mco_l
        mf_to_teff = early_or_late == "early" ? mf_to_teff_e : mf_to_teff_l
        mf_to_mhe  = early_or_late == "early" ? mf_to_mhe_e  : mf_to_mhe_l
        mf_to_xc   = early_or_late == "early" ? mf_to_xc_e   : mf_to_xc_l
        mf_to_mf   = early_or_late == "early" ? mf_to_mf_e   : mf_to_mf_l

        Mco_end       =  interpol_single_Mmax["M_cocore_end"](Mi)*mf_to_mco(Mi, f_acc)
        Mhe_end       =  interpol_single_Mmax["M_hecore_end"](Mi)*mf_to_mhe(Mi, f_acc)
        X_c_endburnHe =  interpol_single_Mmax["X_C"](Mi)         *mf_to_xc(Mi, f_acc) 
        M_end         =  interpol_single_M_pMS["M_end"](Mi)      *mf_to_mf(Mi,f_acc)
        H = M_end - Mhe_end

        if  M_end <= Mhe_end 
            H = 0
            M_end = Mhe_end
        end   


        BSG = mf_to_teff(Mi, f_acc) > 3.9 
        CC = does_it_CC(Mhe=Mhe_end, Mco=Mco_end)

        (!CC) ? (return :WD,  output_run_data(min(M_end,1.4); solver_string=solver_string) ) : nothing  
        

        if isnan(X_c_endburnHe)
            
            @printf("%.3f %.3f %.3f %.3f %.3f\n\n", M, Mi, X_c_endburnHe, Mhe_end, Mco_end)
            throw(ErrorException)

        end

        SN_explosion, exp_param = will_it_explode(;Mhe=Mhe_end, Mco=Mco_end+0.5, 
                                                  Xc= X_c_endburnHe, Mmax=Mi, 
                                                  when="end_BURN_C", criterion=criterion, 
                                                  reference_grid=reference_grid) 
       
        output = output_run_data(M_end, Mhe_end, Mco_end, X_c_endburnHe, caseC ? Inf : -1, exp_param; solver_string=solver_string)

    else
        throw(ErrorException("SN-TYPE error: undefined kind flag - SINGLE-SN"))
    end
    
    Type = nothing
    if SN_explosion
        (H  >= 0.001 && caseC)          ? (Type = :IIn  ) : nothing
        (!caseC && H  >= 1.0 ) && !BSG  ? (Type = :IIP  ) : nothing 
        (!caseC && H  >= 1.0 ) && BSG   ? (Type = :SN87A) : nothing 
        (!caseC && 0.001<= H <1.00)     ? (Type = :IIb  ) : nothing 
        (!caseC && H  <  0.001 )        ? (Type = :Ibc  ) : nothing 
        (H  < 0.001  && caseC) ? (Type = Ibn; throw(ErrorException)) : nothing
        caseC && print("interacting SN!\n")
    else
        Type = :BH
    end   

    return Type, output 

end 



function SN_OUTPUT_manual(M_end, M_he_core_end, M_co_core, Mmax, Mh, Xc, when, criterion;
    II_I_MH_threshold=0.001, deltaM_C = 0, reference_grid = "MW-ZAMS", logTeff = NaN, solver_string="")

    CC = does_it_CC(Mhe=M_he_core_end, Mco=M_co_core)

    (!CC) ? (return :WD, output_run_data(min(M_end,1.4); solver_string=solver_string) ) : nothing 
    SN_explosion, exp_param = will_it_explode(;   Mhe=M_he_core_end, Mco=M_co_core, 
                                                  Xc=Xc, Mmax=Mmax, 
                                                  when=when, criterion=criterion, 
                                                  reference_grid=reference_grid) 
    BSG = logTeff > 3.9

    Type = nothing 

    if SN_explosion
        (Mh  >= 1.0 ) && !BSG ? (Type = :IIP  ) : nothing 
        (Mh  >= 1.0 ) &&  BSG ? (Type = :SN87A) : nothing 
        (0.001<= Mh <1.00)    ? (Type = :IIb  ) : nothing 
        (Mh  <  0.001 )       ? (Type = :Ibc  ) : nothing 

        if deltaM_C > 0
            if Type in [:IIP, :IIb, :SN87A] 
                Type = :IIn
            else 
                Type = :Ibn
            end
        end
    else
        Type = :BH
    end   

    # deltaM_C > 0 && (SN["II-i"]+SN["Ibc-i"] > 0 ? print("yes!\n") : print("NO!\n"))
    return Type, output_run_data(M_end, M_he_core_end, M_co_core, Xc, deltaM_C, exp_param; solver_string=solver_string) 
end

global counter1 = 0
global counter2 = 0

function SN_OUTPUT_firstSN(terminated, f, model, criterion, solverflag1; which="1",
                   II_I_MH_threshold=0.001)

    Mcore_co_end = f[!, "M_co_core_"*which*"_end_"*which][model]
    Mcore_he_end = f[!, "M_he_core_"*which*"_end_"*which][model]
    Mcore_he_hedep = f[!, "M_he_core_"*which*"_end_BURN_He_"*which][model]
    Mh           = f[!, "M_tot_h_"*which*"_end_"*which][model]
    deltaM_C     = which== "1" ? killnan(f.deltaM_C[model]) : -1
    deltaM_A     = f.deltaM_A[model]
    deltaM_B     = f.deltaM_B[model]
    case = "" 
    (f.deltaM_A[model] > 0 || f.M_1_ini_RLOF_A[model] > 0) && (case *="A") 
    (f.deltaM_B[model] > 0 || f.M_1_ini_RLOF_B[model] > 0) && (case *="B") 
    (f.deltaM_C[model] > 0 || f.M_1_ini_RLOF_C[model] > 0) && (case *="C")

    reference_grid = (occursin("A", case) || occursin("B", case)) && !occursin("2", solverflag1) ?  "MW-Stripped" : "MW-ZAMS"

    Mconv_HeBurn = f[!, "Mconv_max_Heburn_"*which][model]
    Rl_hedep = f[!, "Rl_"*which*"_end_BURN_He_"*which][model]
    stripped     = f[!, "age_"*which*"_stripping_"*which][model] > 0
    SN_explosion = nothing
    exp_param = nothing
    M_end = f[!, "M_"*which*"_end_"*which][model]
    logM_t = to_key(f[!, "logM"][model])
    Mmax = which == "1" ? 10 ^f.logM[model] : 10^f.logM[model] * f.q[model]
    Mcore_max = f[!, "M_he_core_"*which*"_end_"*which][model]
    M_env = f[!, "M_1_end_1"][model]-f[!, "M_he_core_1_end_1"][model]
    Xc  = f[!, "Cc_"*which*"_end_BURN_He_"*which][model]

    caseC_CEE = (case == "C") && f.q[model]<0.70
    if caseC_CEE
        Menv = M_end-Mcore_he_end
        deltaM_C += Menv 
        M_end = Mcore_he_end
        Mh = Inf
    end

    Type = nothing

    if occursin("pHeB", terminated) || occursin("CB", terminated) || occursin("Cdep", terminated) || occursin("Hedep", terminated)      
        #IF THE MODEL TERMINATED "SUCCESFULLY" CHECK EXPLODABILITY/WD and then continue

        CC = does_it_CC(Mhe=Mcore_he_end, Mco=Mcore_co_end, Mconvhe=Mconv_HeBurn)

        (!CC) ? (return :WD, output_run_data(min(M_end,1.4); solver_string=solverflag1) ) : nothing  

        reached_endHeburn  = f[!, "age_"*which*"_end_BURN_He_"*which][model] > 0
        reached_iniCburn   = f[!, "age_"*which*"_ini_BURN_C_"*which][model] > 0
        when_end = ( reached_iniCburn ? "_ini_BURN_C_"  : "?")
        (when_end == "?") && (when_end = (reached_endHeburn ? "_end_BURN_He_" : "?"))
        if when_end == "?" 
            @printf("HEY! YOU SAID THE MODEL REACHED A GOOD TERMINATION POINT! AND YET HERE WE FOUND OTHERWISE!")
            @printf("termination flag: %s", terminated)
            print_history(f, model, ["M_1", "M_2", "M_he_core_1", "M_he_core_2", "R_1", "R_2", "Yc_1", "Yc_2"], which)
            throw(ErrorException("ERROR WITH PS20 - issue with declaring the final state"))
        end 

        Mco = f[!, "M_co_core_"*which*when_end*which][model]
        SN_explosion, exp_param = will_it_explode(;   Mhe=Mcore_he_end, Mco=Mco, 
                                                  Xc=Xc, Mmax=Mmax, 
                                                  when=when_end[2:end-1], criterion=criterion, 
                                                  reference_grid=reference_grid, Mconv_HeBurn=Mconv_HeBurn) 
   

    elseif occursin("HeB", terminated) || ( occursin("pMS", terminated))
        if ! (deltaM_A > 0)  &&  !(deltaM_B >0)
            logM = f[!, "logM"][model]
            q = f[!, "q"][model]
            RL = f[!, "Rl_"*which*"_end_"*which][model]
            R_max   =  f_Rmax_CB(10 .^ logM * (which == "1" ? 1 : q))
            R_caseB = f_Rmax_pMS(10 .^ logM * (which == "1" ? 1 : q))

            if RL >= R_max #is a single star
                solverflag1 *= "no RLOF (extrap) - SINGLE_SN"
                return SINGLE_SN(Mmax, criterion, "single"; solver_string=solverflag1)
            elseif RL >= R_caseB  #would undergo RLOF - not case B but case C -> interacting SN

                # output = SINGLE_SN(SN, Mmax, criterion, "single"; solver_string=solverflag1)
                Mi = 10^logM
                M_end = interpol_single_Mmax["M_end"](Mi)
                Mcore_he_end = interpol_single_Mmax["M_hecore_end"](Mi)
                Mcore_co_end = interpol_single_Mmax["M_cocore_end"](Mi)
                Xc = interpol_single_Mmax["X_C"](Mi)
                Menv = M_end-Mcore_he_end
                if q < 0.70
                    deltaM_C = Menv 
                    M_end = Mcore_he_end
                    Mh = Inf
                end
                
                Xc = interpol_single_Mmax["Xc"](10^logM)
                solverflag1 *= "Case C (extrap) - Manual + inter"

                Type, output = SN_OUTPUT_manual(M_end, Mcore_he_end, Mcore_co_end, Mi, Mh, Xc, "end_BURN_C", criterion, reference_grid = "MW-Stripped"; solver_string=solverflag1)

                if Type in [:IIP, :IIb]
                    Type=:IIn
                end  

                return Type, output
            else #should have undergone case B after all! 
                if !occursin("pMS", terminated) 
                    print(terminated, "\n")
                    print_history(f, model, ["M_1", "M_he_core_1", "M_2", "M_he_core_2",  "R_1", "Rl_1", "Rl_2", "R_2"], which)
                    throw(ErrorException("a HeB that should have undergone Case B yet it didn't?"))
                end 

                solverflag1 *= "Case B (extrap) - SINGLE_SN (stripped)"

                Type, output = SINGLE_SN(Mcore_max, criterion, "stripped"; solver_string=solverflag1)
                return Type, output
            end
        else
            if f[!, "Yc_1_end_1"][model]<0.25 
                Mco = interpol_HeS_Mmax["M_cocore_end"](Mcore_max)
                Xc = interpol_HeS_Mmax["X_C"](Mcore_max)

                solverflag1 *= "Case B (extrap) - SN_OUTPUT_manual (stripped)"

                Type, output = SN_OUTPUT_manual(M_end, Mcore_max, Mco, Mmax, M_env, Xc, "end_BURN_He", criterion, reference_grid = "MW-Stripped"; solver_string=solverflag1)
                return Type, output
            else 
                #OBSERVED BY EYE! MOST OF THESE MODELS ARE VERY MASSIVE HESTARS
                #WITH THIN-ish H-RICH ENVELOPES 

                solverflag1 *= "outliers - SINGLE_SN (stripped)"

                Type, output = SINGLE_SN(Mcore_max, criterion, "stripped"; solver_string=solverflag1)
                #print_history(f, model, ["M_1", "M_he_core_1", "Yc_1", "M_2", "M_he_core_2", "Yc_2", "R_1", "Rl_1", "Rl_2", "R_2"], which)
                #print(SN)
                return Type, output
            end
        end
        throw(ErrorException("SHOULD HAVE RETURNED ALREADY!"))

    else
        solverflag1 *= " early termination - SN_OUTPUT_firstSN->SingleSN "

        Type, output = SINGLE_SN(10^f.logM[model], criterion, "single")
        return Type, output 


        @printf("model: %5d,   Star = %s\n", model ,which)
        @printf("Terminated = %s\n", terminated)
        @printf("ERROR! Terminated too early.\n")
        @printf("Termination string: %s\n ", f.summary[model])
        c = 0
        for col in names(f)
            c+=1
            str = strip(col)
            ! (col in ["summary", "STATUS"]) && @printf("%30s  =  %10.4e %7s", str, f[!, str][model], "")
            (c%3==0) ? @printf("\n") : nothing 
        end

        throw(ErrorException("SN-TYPE error: undefined outcome - early termination"))
    end

    solverflag1 *= " nominal"

    if !SN_explosion
        Type = :BH 
        return Type, output_run_data(M_end; solver_string=solverflag1) 
    end


    
    Hrich =  (Mh >= II_I_MH_threshold) 
    narrow_lines =  deltaM_C > 0

    if Hrich
        Type =  (Mh >= 1)  ? :IIP : :IIb
        narrow_lines && (Type = :IIn)
    else
        dM=0
        if  predictor_caseBC(f, model; which=which)
            dt = retrieve_tpreCC(Mconv_HeBurn, Rl_hedep)/1e3
            # dM = 0.789*(1-exp(-(dt/5.22)^1.31))
            dM = 0.769*(1-exp(-(dt/5.49)^1.84))
            if dM > 0
                deltaM_C = dM
                narrow_lines = deltaM_C > 0 
                M_end = Mcore_he_hedep-dM 
                M_end -= 0.1212 #extra 0.15 Msun to account for the wind that is not included
                Mcore_he_end = M_end
            end
        end 
        Type = narrow_lines ? :Ibn : :Ibc
    end

    return Type, output_run_data(M_end, Mcore_he_end, Mcore_co_end, Xc, deltaM_C, exp_param; solver_string=solverflag1) 
end



function SN_OUTPUT_secondSN(terminated, f, model, criterion, kick, which, firstSN, endvals_1stSN, pre_SN_orbit, post_SN_orbit, solverflag2)

    entry(what, of_which, when, whose) = f[!, what * "_" * of_which * "_" * when * "_" * whose][model]
    this_star = which 
    the_star_which_went_SN_first = ( which == "1") ? "2" : "1"
    logM = f.logM[model]
    q = f.q[model]
    M0 = entry("M", this_star, "end", the_star_which_went_SN_first)
    isnan(M0) && @printf("WARNING! M0=NaN! model=%d", model)
    Mmax = max(M0, this_star == "1" ? 10^logM : 10^logM * q )
    extraMT = ""

    Type = nothing 

    M_at(when) = f[!, "M_" * this_star * "_" * when * "_" * this_star][model]
    Mh_at(when) = f[!, "M_tot_h_" * this_star * "_" * when * "_" * this_star][model]
    R_at(when) = f[!, "R_" * this_star * "_" * when * "_" * this_star][model]
    t_at(when) = f[!, "age_" * this_star * "_" * when * "_" * this_star][model]
    t_SN1 = entry("age", the_star_which_went_SN_first, "end", the_star_which_went_SN_first)

    if (post_SN_orbit.orbit==:unbound) && !(firstSN == :BH  || firstSN == :WD)
        solverflag2 *= "Unbound-"
        if occursin("pHeB", terminated) || occursin("CB", terminated) || occursin("Cdep", terminated)
            solverflag2 *= "SN_OUTPUT_MANUAL"  

            occursin("pHeB", terminated) && (when_end = "end_BURN_He")
            occursin("CB",   terminated) && (when_end = "ini_BURN_C")
            occursin("Cdep", terminated) && (when_end = "end_BURN_C")
            M_end = f[!, "M_" * this_star * "_"*when_end*"_" * this_star][model]
            M_he_core_end = f[!, "M_he_core_" * this_star * "_"*when_end*"_" * this_star][model]
            M_co_core = f[!, "M_co_core_" * this_star * "_"*when_end*"_" * this_star][model]
            Menv = M_end-M_he_core_end
            Xc = f[!, "Cc_" * this_star * "_end_BURN_He_" * this_star][model]
            Type, output =  SN_OUTPUT_manual(M_end, M_he_core_end, M_co_core, Mmax, Menv, Xc, when_end, criterion;
                                        reference_grid = "MW-ZAMS", solver_string=solverflag2)
            return Type, output, extraMT
        else 
            solverflag2 *= "SINGLE_SN"  
            Type, output =  SINGLE_SN(Mmax, criterion, "single"; solver_string=solverflag2)
            return Type, output, extraMT
        end 
        
    elseif post_SN_orbit.orbit == :bound || (firstSN == :BH  || firstSN == :WD)
        Rl = nothing
        solverflag2 *= "Bound-"



        #determine post-SN RL
        if firstSN == :BH
            solverflag2 *= "RL(BH)-"
            #nothing happens - new RL is the same as before the explosion
            Rl = entry("Rl", this_star, "end", the_star_which_went_SN_first) 
        elseif firstSN == :WD
            solverflag2 *= "RL(WD)-"
            Rl = eval_RL(entry("M", this_star, "end", the_star_which_went_SN_first), endvals_1stSN.M_end,       "a", pre_SN_orbit.a)
        else
            solverflag2 *= "RL(NS)-"
            Rl = eval_RL(entry("M", this_star, "end", the_star_which_went_SN_first), endvals_1stSN.M_remnant_g, "a", post_SN_orbit.a_peri/Rsun)
        end

        #determine post-SN R history
        logM = log10(entry("M", this_star,  "end", the_star_which_went_SN_first) )
        logP = log10(f[!, "period_end_"*the_star_which_went_SN_first][model])

        #CHECK TAMS -> CASE A
        R = occursin("MS",terminated) ?  f_Rmax_MS(Mmax) : R_at("end_BURN_H")
        if (R >= Rl)  #case A RLOF -> likely TZO -> no SN 
            extraMT = "A"
            solverflag2 *= "dissolved"  
            output = output_run_data(0; solver_string=solverflag2) 
            return :X, output, extraMT
        end 

        #CHECK ini_BURN_He -> CASE B
        R = ( occursin("pMS",terminated) || occursin("MS", terminated)) ? f_Rmax_pMS(Mmax) :  R_at("ini_BURN_He")
        if (R >= Rl) #case B RLOF
            #assume CE ejection - short binary 
            extraMT = "eB"
            Mh = 0 
            Mcore = nothing 
            if occursin("MS",terminated)
                Mcore = interpol_single_Mmax["M_hecore_pMS"](Mmax) 
            elseif occursin("pMS",terminated)
                Mcore = entry("M_he_core", this_star, "end", this_star)
            else 
                Mcore = entry("M_he_core", this_star, "ini_BURN_He", this_star)
            end

            isnan(Mcore) && throw(ErrorException)
            Mhecore = interpol_HeS_Mmax["M_hecore_end"](Mcore)
            Mco =interpol_HeS_Mmax["M_cocore_end"](Mcore)
            Xc = interpol_HeS_Mmax["X_C"](Mcore)        
           
            solverflag2 *= @sprintf("caseB (extrap: R vs Rl : %.2e-%.2e) - SN_OUTPUT_manual - ", R, Rl)  

            Type, output = SN_OUTPUT_manual(Mhecore, Mhecore, Mco, Mmax, 0., Xc, "ini_BURN_He", criterion, reference_grid = "MW-Stripped"; solver_string=solverflag2)                
            Type==:Ibn && throw(ErrorException)

            return Type, output, extraMT
        end 

        #CHECK if it would RLOF during core HeB -> "late" CASE B
        R = max(R_at("end"), f_Rmax_HeB(Mmax))
        if  ( occursin("pMS",terminated) || occursin("MS", terminated) || occursin("HeB",terminated)) &&  R >= Rl 
            #assume CE ejection - short binary 
            extraMT = "lB"

            Mh = ( occursin("pMS",terminated) || occursin("HeB",terminated)) ? Mh_at("end") : 1e99
            Mhecore   = Mh > 1. ? interpol_single_Mmax["M_hecore_end"](Mmax) : interpol_single_Mmax["M_hecore_pMS"](Mmax)
            Mco = interpol_HeS_Mhedep["M_cocore_end"](Mhecore)
            Xc = interpol_HeS_Mmax["X_C"](Mhecore)
            deltaM_C = 0
            solverflag2 *= "caseB - SN_OUTPUT_manual"  
            Type, output = SN_OUTPUT_manual(Mhecore, Mhecore, Mco, Mmax, 0., Xc, "end_BURN_He", criterion, reference_grid = "MW-ZAMS", deltaM_C=deltaM_C; solver_string=solverflag2)          
            return Type, output, extraMT

        end
        
        #check if it would RLOF after core HeB
        R = max(R_at("end"), f_Rmax_CB(Mmax))

        if  ( occursin("pMS",terminated) || occursin("MS", terminated) || ( occursin("HeB",terminated) && !occursin("pHeB",terminated)) ) &&  R >= Rl 
            #assume CE ejection - short binary 

            extraMT = "C"

            Mh = ( occursin("pMS",terminated) || occursin("HeB",terminated)) ? Mh_at("end") : Inf
            Mhecore   = Mh > 1. ? interpol_single_Mmax["M_hecore_end"](Mmax) : interpol_single_Mmax["M_hecore_pMS"](Mmax)
            Mco = interpol_HeS_Mhedep["M_cocore_end"](Mhecore)
            Xc = interpol_HeS_Mmax["X_C"](Mhecore)
            deltaM_C = Inf
            solverflag2 *= "caseC (extrap) - SN_OUTPUT_manual"  
            Type, output = SN_OUTPUT_manual(Mhecore, Mhecore, Mco, Mmax, 0., Xc, "end_BURN_He", criterion, reference_grid = "MW-ZAMS", deltaM_C=deltaM_C; solver_string=solverflag2)          
            Type==:Ibn && @printf("termination %s , Mh %5.3f  Mhecore %5.3f Mcocore %5.3f \n", terminated, Mh, Mhecore, Mco, )
            Type==:Ibn && @printf("Ibc-i! %5s -> Rmax vs Rl = %.1f - %.1f sec pMs/Ms/HeB\n", model, R, Rl)
            Type==:Ibn && print_history(f, model, ["M_1", "M_2", "M_he_core_1", "M_he_core_2", "R_1", "R_2", "Yc_1", "Yc_2",  "Rl_1", "Rl_2",], which)
            
            return Type, output, extraMT

        elseif (R_at("end") >= Rl)  
            extraMT = "C"
            #assume CE ejection - short binary 
            M = M_at("end")
            Mhecore = entry("M_he_core", this_star, "end", this_star)
            Mcocore = entry("M_co_core", this_star, "end_BURN_He", this_star)
            Xc = entry("Cc", this_star, "end_BURN_He", this_star)
            Menv = M-Mhecore
            deltaM_C = Inf
            output = nothing 
            solverflag2 *= "caseC - SN_OUTPUT_manual"  
            try 
                Type, output = SN_OUTPUT_manual(M, Mhecore, Mcocore, Mmax, Menv, Xc, "ini_BURN_C", criterion, reference_grid = "MW-ZAMS", deltaM_C=deltaM_C; solver_string=solverflag2)                
            catch 
                @printf("termination of star %s : %s\n", which, terminated)
                print_history(f, model, ["M_1", "M_2", "M_he_core_1", "M_he_core_2", "M_co_core_1", "M_co_core_2", "R_1", "R_2", "Yc_1", "Yc_2", "Rl_1", "Rl_2"], which)
                throw(ErrorException("SN EXPLOSION FAILED!"))
            end


            # SN["Ibc-i"]>0 && @printf("Ibc-i! %5s -> Rmax vs Rl = %.1f - %.1f sec end\n", model, R, Rl)

            return Type, output, extraMT
        end 
        
        if occursin("pHeB", terminated) || occursin("CB", terminated) || occursin("Cdep", terminated)  

            solverflag2 *= "no RLOF - SN_OUTPUT_manual"  
            occursin("pHeB", terminated) && (when_end = "end_BURN_He")
            occursin("CB",   terminated) && (when_end = "ini_BURN_C")
            occursin("Cdep", terminated) && (when_end = "end_BURN_C")
            M_end = f[!, "M_" * this_star * "_"*when_end*"_" * this_star][model]
            M_he_core_end = f[!, "M_he_core_" * this_star * "_"*when_end*"_" * this_star][model]
            M_co_core = f[!, "M_co_core_" * this_star * "_"*when_end*"_" * this_star][model]
            Menv = M_end-M_he_core_end
            Xc = f[!, "Cc_" * this_star * "_end_BURN_He_" * this_star][model]
            
           Type, output =  SN_OUTPUT_manual(M_end, M_he_core_end, M_co_core, Mmax, Menv, Xc, 
                                        when_end, criterion; 
                                        reference_grid = "MW-ZAMS", solver_string=solverflag2)
            return Type, output, extraMT
        else 
            solverflag2 *= "no RLOF - SINGLE_SN"  
            Type, output = SINGLE_SN(Mmax, criterion, "single"; solver_string=solverflag2)

            return Type, output, extraMT
        end 


        @printf("OH NO! THE SECOND SN HERE WAS NOT CAPTURED CORRECTLY!\n")
        @printf("termination flag: %s\n", terminated)
        @printf("Rl: %7.2f\n", Rl)
        @printf("Mmax: %7.2f\n", Mmax)
        @printf("Rmax(MS): %7.2f\n", f_Rmax_MS(Mmax))
        @printf("Rmax(pMS): %7.2f\n", f_Rmax_pHeB(Mmax))
        @printf("Rmax(end): %7.2f\n", f_Rmax_CB(Mmax)  )
        print_history(f, model, ["M_"*this_star, "R_"*this_star], this_star)
        throw(ErrorException("THE SECOND SN HERE WAS NOT CAPTURED CORRECTLY"))

    elseif (post_SN_orbit.orbit==:direct_collision) 
        solverflag2 *= "direct_collision-"
        solverflag2 *= "dissolved"  
        output = output_run_data(0; solver_string=solverflag2) 
        return :X, output, extraMT
        
    else 
        print(kick, post_SN_orbit)
        throw(ErrorException("unrecognized kick!"))
    end

end


function is_merger(f, model, criterion, terminated1, terminated2)
    merger = false 
    M1 = NaN
    M1core = NaN
    M2 = NaN 
    M2core = NaN
    M1h = NaN 
    M2h = NaN
    t_merger = NaN 
    Xc = NaN
    label = "no merger"
    case = nothing 
    early_or_late = ""
    logT_threshold_caseB = 3.6 
    Xc_threshold_caseA = 0.15
    logM = f.logM[model]
    logP = f.logP[model]
    q = f.q[model]


    if terminated1 == "not run: ZAMS merger" && terminated2 == "not run: ZAMS merger"
        return Dict("merger?"=>true, 
            "M1"=> 10. ^logM, "M2"=> 10. ^logM * q , 
            "M1core"=>0., "M2core"=>0.,
            "M1h"=>10. ^logM, "M2h"=>10. ^logM * q,
            "Xc" => 0.70,
            "t_merger"=>NaN, "case"=>"A-mer", "label"=>"case A not run (ZAMS)", "early_or_late"=>"early",
            "RLOFcase" => "X", "CaseC_CEE"=>false)
        end 
    
    if (! isnan(f.deltaM_A[model])) && isnan(f.age_1_end_BURN_H_1[model])
        merger = true 
        label = "case A unfinished"
        case = "A-mer-unf"
        early_or_late = (f.Xc_1_ini_RLOF_A[model] <= Xc_threshold_caseA) ? "late" : "early"

    elseif (! isnan(f.deltaM_B[model])) && (isnan(f.age_1_end_RLOF_B[model]) || f.age_1_end_1[model]  < f.age_1_end_RLOF_B[model]+10 )
        merger = true 
        label = "case B unfinished"
        case = "B-mer-unf"
        early_or_late = (f.logT_1_ini_RLOF_B[model] <= logT_threshold_caseB) ? "late" : "early"
    elseif  (f.R_2_end_1[model] >= 0.95 * f.Rl_2_end_1[model]) || f.STATUS[model] == "IMT"
        #@printf("%4d: R2 vs Rl2 @end1: %6.1f - %6.1f\n", model, f.R_2_end_1[model], f.Rl_2_end_1[model])
        merger = true 
        case = "IMT"
        label = "Inverse Mass Transfer"
        if isnan(f.age_1_end_BURN_H_1[model])
            case*="-A(1)"
            early_or_late = (f.Xc_1_end_1[model] <= Xc_threshold_caseA) ? "late" : "early"
        elseif isnan(f.age_1_end_BURN_He_1[model])
            case*="-B(1)"
            early_or_late = (f.logT_1_end_1[model] <= logT_threshold_caseB) ? "late" : "early"
        elseif isnan(f.age_1_ini_BURN_C_1[model])
            case*="-C(1-early)"
            # (f.M_co_core_1_end_1[model] > 1.4) && print(to_key(logM)*"/"*to_key(logP)*"/"*to_key(q)*" ", model, " ", case,  " This is a candidate for an interacting SN!", terminated1, " ", terminated2, "\n")
        elseif isnan(f.age_1_end_BURN_C_1[model])
            case*="-C(1-late)"
            # (f.M_co_core_1_end_1[model] > 1.4) && print(to_key(logM)*"/"*to_key(logP)*"/"*to_key(q)*" ",model, " ", case,  " This is a candidate for an interacting SN!", terminated1, " ", terminated2, "\n")
        else 
            throw(ErrorException("IMT - stage unclear!"))
        end
    elseif (f.R_1_end_1[model] >= 0.95 * f.Rl_1_end_1[model]) 
        if isnan(f.age_1_end_BURN_H_1[model])
            merger = true
            case="A-mer"
           #@printf("signalled merger [A](%s-%s)\n", terminated1, terminated2,)
            label = "Terminated while close to (or during) RLOF-"*case[1]
            early_or_late = (f.Xc_1_end_1[model] <= Xc_threshold_caseA) ? "late" : "early"

        elseif occursin("pMS", terminated1)  
             case="B-mer"
             merger = true
             #@printf("signalled merger [B](%s-%s)\n", terminated1, terminated2,)
             label = "Terminated while close to (or during) RLOF-"*case[1]
             early_or_late = (f.logT_1_end_1[model] <= logT_threshold_caseB) ? "late" : "early"

        else 
            if occursin("caseC", terminated1) || occursin("pHeB", terminated1) || occursin("CB", terminated1) || occursin("Cdep", terminated1)
                #the model likely goes through case C - ignore for merger!
            else
                @printf("%5s %10s %10s - (cor %.2f/env %5.2f) %5.0f - %5.0f | %5.0f - %5.0f\n", model, terminated1,  terminated2, f.M_he_core_1_end_1[model],  f.M_1_end_1[model]-f.M_he_core_1_end_1[model], f.R_1_end_1[model], f.Rl_1_end_1[model],  f.R_2_end_1[model], f.Rl_2_end_1[model])
                throw(ErrorException("this model is close to RLOF but crashed -> merger\nhowever the end is unclear"))
            end
        end
        # label = "Terminated while close to (or during) RLOF-"*case[1]
    end
    

    #if !merger  #check summary column!
    #     summary_string = f[!, "summary"][model]
    #     merger_flags = [ 
    #         "Terminated_due_to_mass_transfer_rate_reaching_maximum", 
    #         "ZAMS_L2_overflow",
    #     ]
    #     for flag in merger_flags
    #         if occursin(flag, summary_string)
    #             merger = true 
    #             break
    #         end
    #     end
    # end

    if merger 
        M1 = f.M_1_end_1[model]
        M2 = f.M_2_end_2[model]
        M1core = f.M_he_core_1_end_1[model]
        M2core = f.M_he_core_2_end_2[model]
        t_merger = f.age_1_end_1[model]
        M1h = f.M_tot_h_1_end_1[model] 
        M2h = f.M_tot_h_2_end_2[model]
        Xc = f.Xc_1_end_1[model]
        early_or_late = f.logT_1_end_1[model] <= 3.6 ? "late" : "early" 
    end

    for rlof in ["A", "B"]
        if criterion in ["ERK", "PA_IV", "PAULI", "PABLO"]
            extra_string ="_unstable_RLOF_"*rlof*"_"*criterion
            MERGER = f[!, "age_1"*extra_string][model]>0
            if MERGER

                merger = true 
                M1 = f[!, "M_1"*extra_string][model]
                M2 = f[!, "M_2"*extra_string][model]
                M1core = f[!, "M_he_core_1"*extra_string][model]
                M2core = f[!, "M_he_core_2"*extra_string][model]
                M1h = f[!, "M_tot_h_1"*extra_string][model] 
                M2h = f[!, "M_tot_h_2"*extra_string][model]
                t_merger = f[!, "age_1"*extra_string][model]
                Xc = f[!, "Xc_1"*extra_string][model]
                label = "Unstable Case "*rlof*" per "*criterion
                case = rlof*"-merC"
                if rlof == "A" 
                    early_or_late = f[!, "Xc_1_ini_RLOF_"*rlof][model] <= Xc_threshold_caseA ? "late" : "early" 
                elseif rlof == "B"
                    early_or_late = f[!, "logT_1_ini_RLOF_"*rlof][model] <= logT_threshold_caseB ? "late" : "early" 
                else
                    early_or_late = "" 
                end

                break
            end
        end
    end 

    rlof = ""
    (f.deltaM_A[model] > 0 || f.M_1_ini_RLOF_A[model] > 0) && (rlof *="A") 
    (f.deltaM_B[model] > 0 || f.M_1_ini_RLOF_B[model] > 0) && (rlof *="B") 
    (f.deltaM_C[model] > 0 || f.M_1_ini_RLOF_C[model] > 0) && (rlof *="C") 

    if criterion == "ALL" && rlof != "" && rlof != "C"
        extra_string ="_ini_RLOF_"*rlof[1]
        merger = true 
        M1 = f[!, "M_1"*extra_string][model]
        M2 = f[!, "M_2"*extra_string][model]
        M1core = f[!, "M_he_core_1"*extra_string][model]
        M2core = f[!, "M_he_core_2"*extra_string][model]
        M1h = f[!, "M_tot_h_1"*extra_string][model] 
        M2h = f[!, "M_tot_h_2"*extra_string][model]
        t_merger = f[!, "age_1"*extra_string][model]
        Xc = f[!, "Xc_1"*extra_string][model]
        label = "Unstable Case "*rlof*" per "*criterion
        case = rlof[1]*"-mer-all"
        if occursin("A",rlof) 
            early_or_late = f[!, "Xc_1_ini_RLOF_"*rlof[1]][model] <= Xc_threshold_caseA ? "late" : "early" 
        elseif occursin("B",rlof) 
            early_or_late = f[!, "logT_1_ini_RLOF_"*rlof[1]][model] <= logT_threshold_caseB ? "late" : "early" 
        else
            early_or_late = "" 
        end

    end


    #HARDCODED! 
    #UNSTABLE CASE C IS SUCH ONLY FOR CASE C SYSTEMS WITH q<0.7.
    CaseC_CEE = (rlof == "C") && (f.q[model] < 0.70 || criterion == "ALL")

    return Dict("merger?"=>merger, 
                "M1"=>M1, "M2"=>M2, 
                "M1core"=>M1core, "M2core"=>M2core,
                "M1h"=>M1h, "M2h"=>M2h,
                "Xc" => Xc,
                "t_merger"=>t_merger, "case"=>case, "label"=>label, "early_or_late"=>early_or_late, "RLOFcase" => rlof, "CaseC_CEE"=>CaseC_CEE)
end

function not_run(f, model)
    return f.STATUS[model] == "n.a."
    # return isnan(f.age_1_end_BURN_H_1[model]) && isnan(f.age_2_end_BURN_H_2[model]) && isnan(f.deltaM_A[model])
end

function termination_stage(f, model, which)
       
    terminated = ""
    entry(what, when) = f[!, what*"_"*which*"_"*when*"_"*which][model]
    close_to_RLOF = entry("R", "end") > 0.95 * entry("Rl", "end")
    if isnan(entry("age", "end_BURN_H")) && entry("Xc", "end") > 1e-6
        terminated *= "MS"
        if f[!, "age_"*which*"_ini_RLOF_A"][model] > 0 && isnan(f[!, "age_"*which*"_end_RLOF_A"][model]) || close_to_RLOF
            terminated *= "-caseA"
        end
    elseif isnan(entry("age", "ini_BURN_He"))
        terminated *= "pMS"
        if f[!, "age_"*which*"_ini_RLOF_B"][model] > 0 && isnan(f[!, "age_"*which*"_end_RLOF_B"][model])
            terminated *= "-caseB"
        end
    elseif isnan(entry("age", "end_BURN_He")) && entry("Yc", "end") > 1e-6
        terminated *= "HeB"
        if f[!, "age_"*which*"_ini_RLOF_B"][model] > 0 && isnan(f[!, "age_"*which*"_end_RLOF_B"][model])
            terminated *= "-caseB"
        end
    elseif isnan(entry("age", "end_BURN_He")) && entry("Yc", "end") <= 1e-6
        terminated *= "Hedep"
        if f[!, "age_"*which*"_ini_RLOF_B"][model] > 0 && isnan(f[!, "age_"*which*"_end_RLOF_B"][model])
            terminated *= "-caseB"
        end
    elseif isnan(entry("age", "ini_BURN_C")) 
        terminated *= "pHeB"
        if f[!, "age_"*which*"_ini_RLOF_C"][model] > 0 || f.mdot_RLOF_max_RLOF_C[model] > 1e-5
            terminated *= "-caseC"
        end
    elseif isnan(entry("age", "end_BURN_C"))
        terminated *= "CB"
        if f[!, "age_"*which*"_ini_RLOF_C"][model] > 0  || f.mdot_RLOF_max_RLOF_C[model] > 1e-5
            terminated *= "-caseC"
        end
    else
        terminated = "Cdep"
        if f[!, "age_"*which*"_ini_RLOF_C"][model] > 0  || f.mdot_RLOF_max_RLOF_C[model] > 1e-5
            terminated *= "-caseC"
        end
    end


    return terminated 
end

function stage_when_other_terminated(f, model, which)
    star_output = which
    star_reference = (which=="1") ? "2" : "1"
    terminated = ""
    age(when) = f[!, "age_"*star_output*"_"*when*"_"*star_output][model]
    age_ref = f[!,   "age_"*star_reference*"_end_"*star_reference][model] 
    if     isnan(age("end_BURN_H"))   || age_ref < age("end_BURN_H")
        terminated *= "MS"
    elseif isnan(age("ini_BURN_He"))  || age_ref < age("ini_BURN_He")
        terminated *= "pMS"
    elseif isnan(age("end_BURN_He"))  || age_ref < age("end_BURN_He")
        terminated *= "HeB"
    elseif isnan(age("ini_BURN_C"))   || age_ref < age("ini_BURN_C")
        terminated *= "pHeB"
    elseif isnan(age("end_BURN_C"))   || age_ref < age("end_BURN_C")
        terminated *= "CB"
    else
        terminated = "Cdep"
    end


    return terminated 
end



evaluate_period(M1,M2,a) = 2*pi*sqrt(a^3 / (G*(M1+M2)))
EK(M,V) = 0.5*M*norm(V)^2
EG(M1,M2,R) = -G*M1*M2/norm(R)
is_orbit_unbound(Mold,Mnew) = Mnew <= 0.5*Mold

function derive_a_massloss_smooth(a_pre, M1, M2, DeltaM; dm=1e-4)
    a_old = a_pre
    a_new = a_pre

    if DeltaM > 1e-4 
        for i in range(1, round(Int, DeltaM/dm))
            DM = dm * i
            q = (M1-DM)/M2
            a_new = a_old *(1-2*dm/M1 * (1 - (2+q)/(2*q*(1+q))))
            a_old = a_new
        end 
    end

    m1 = M1-DeltaM
    a_omega_orb  = sqrt(G*(m1+M2)/a_new * Msun/Rsun)
    v2 = a_omega_orb * m1/(m1+M2) / kms
    v1 = a_omega_orb * M2/(m1+M2) / kms

    return (a=a_new, v1=v1, v2=v2)
end 




function get_multiplicity(method; mi=0.01, mf=100., custom_value = NaN)
    # f_mult =CSV.read(data_dir*"Offner2023_multiplicity.dat",
    #                 header=1, DataFrame,delim=' ',ignorerepeated=true)
    # M_mult = 10 .^ ( (log10.(f_mult.M_low) + log10.(f_mult.M_high)) ./ 2)
    # #sM_mult_low, sM_mult_high = M_mult - f_mult.M_low  , f_mult.M_high - M_mult
    # MF = 0.01 * f_mult.MF
    # f_intp = interpolator(log10.(M_mult), MF, extrapolation=:flat)
    # if method == :interpolate
    #     f1(x) =  f_intp( log10.(x) )
    #     return f1
    # elseif method == :avg
    #     logms = range(log10(mi), log10(mf), step = 0.01)
    #     f2(x) = mean(f_intp(logms))
    #     return f2
    if method == :custom
        f5(x) = custom_value
        return f5
    else
        error("Invalid method: $method")
    end


end 

function load_models(filename, MER_CRIT, MER_dM, MER_EOL, EXP_CRIT)

    logMs = range(0.70, 2.00, step=0.05)
    qs = range(0.05, 0.95, step=0.05)
    logPs = range(0.00, 3.75, step=0.05)
    MODELS = SortedDict()
    for logM in logMs 
        MODELS[to_key(logM)] = SortedDict()
        for logP in logPs 
            MODELS[to_key(logM)][to_key(logP)] = SortedDict()
            for q in qs
                MODELS[to_key(logM)][to_key(logP)][to_key(q)] = Dict()
            end
        end
    end
    
    fp = CSV.read(filename, header=1, DataFrame,delim=' ',ignorerepeated=true)   
    
    for i in range(1,length(fp.logM))
        (fp.MER_CRIT[i] != MER_CRIT) ?  continue : nothing
        (fp.MER_dM[i]   != MER_dM  ) ?  continue : nothing
        (fp.MER_EOL[i]  != MER_EOL ) ?  continue : nothing
        (fp.EXP_CRIT[i] != EXP_CRIT) ?  continue : nothing 
        logM = fp.logM[i]
        logP = fp.logP[i]
        q = fp.q[i]
        #print(logM, typeof(logP), "\n")
        #print(to_key(logM), "\n")
        logM_t = to_key(logM)
        logP_t = to_key(logP)
        q_t = to_key(q)

        MODELS[logM_t][logP_t][q_t]["mod"] = fp.mod[i]
        MODELS[logM_t][logP_t][q_t]["M1max"] = fp[!, "M1_max"][i]
        MODELS[logM_t][logP_t][q_t]["M2max"] = fp[!, "M2_max"][i]
        MODELS[logM_t][logP_t][q_t]["M1f"] = fp[!, "M1_end"][i]
        MODELS[logM_t][logP_t][q_t]["M2f"] = fp[!, "M2_end"][i]
        MODELS[logM_t][logP_t][q_t]["t_SN1"] = fp.t_SN1[i]
        MODELS[logM_t][logP_t][q_t]["t_SN2"] = fp.t_SN2[i]
        MODELS[logM_t][logP_t][q_t]["first"] = fp.first[i]
        MODELS[logM_t][logP_t][q_t]["second"] = fp.second[i]
        MODELS[logM_t][logP_t][q_t]["casus"] = fp.casus[i]
        MODELS[logM_t][logP_t][q_t]["SN1"] = SortedDict()
        MODELS[logM_t][logP_t][q_t]["SN2"] = SortedDict()
        for out in outcomes
            key1 = @sprintf("%s-1", out)
            key2 = @sprintf("%s-2", out)
            MODELS[logM_t][logP_t][q_t]["SN1"][out] = fp[!, key1][i]
            MODELS[logM_t][logP_t][q_t]["SN2"][out] = fp[!, key2][i]
        end
    end 

    return MODELS 
end 

outcomes = ["WD", "IIP", "IIP-B", "IIb", "Ibc", "II-i", "Ibc-i", "BH", "X"]
SN_outcomes = copy(outcomes)
filter!(!=("WD"),SN_outcomes)
filter!(!=("X"), SN_outcomes)
filter!(!=("BH"),SN_outcomes)
progenitors = ["1", "2", "M", "S"]


function unstable_mass_transfer_analysis(f, model, criterion, massloss; val = 1)

    a_pre = NaN
    Mcore = NaN
    Rcore = NaN
    M1 =    NaN
    M2 =    NaN
    R2 =    NaN
    Ebind = NaN

    function f_q_RL(q)
        qq = q .^ (-1/3)
        return 0.49 .* qq .* qq ./ (0.6 .* qq .* qq .+ log.(1 .+ qq))
    end 

    if ! (criterion in ["X", "ANY", "ALL"])
        a_pre = f[!, "a_unstable_RLOF_B_"*criterion][model]
        Mcore = f[!, "M_he_core_1_unstable_RLOF_B_"*criterion][model]
        Rcore = f[!, "CE_core_R_1_unstable_RLOF_B_"*criterion][model]
        M1 = f[!, "M_1_unstable_RLOF_B_"*criterion][model]
        M2 = f[!, "M_2_unstable_RLOF_B_"*criterion][model]
        R2 = f[!, "R_2_unstable_RLOF_B_"*criterion][model]
        Ebind = f[!, "Ebind_b_1_unstable_RLOF_B_"*criterion][model]
    end
    if isnan(a_pre)
        a_pre = f.a_end_1[model]
        Mcore = f.M_he_core_1_end_1[model]
        Rcore = f.CE_core_R_1_end_1[model]
        M1 =    f.M_1_end_1[model]
        M2 =    f.M_2_end_1[model]
        R2 =    f.R_2_end_1[model]
        Ebind = f.Ebind_b_1_end_1[model]
    end
    if Mcore <0.001
        return 0
    end 

    if massloss == "Energy"
        a_min  =  max( R2/f_q_RL(Mcore/M2), Rcore/f_q_RL(M2/Mcore) )
        dE_orb = abs(G*Msun*Msun/(2*Rsun) * (Mcore*M2/a_min - M1*M2/a_pre))
        if (isnan(f_q_RL(Mcore/M2))) 
            @printf("mod, Mc, M2, f_q : %5d %5.3f %5.3f %5.3f", model, Mcore, M2, f_q_RL(Mcore/M2))
            throw(ErrorException("mmhmhm"))
        end
        dM_rel_donor = abs(dE_orb/Ebind)

        return dM_rel_donor
    elseif massloss == "FIX"
        return val
    else
        return 0
    end
end 


function print_history(f, model, what, which)
    stages = ["ini_BURN_H_1",
              "BURN_H_50_1",
              "BURN_H_25_1",
              "end_BURN_H_1",
              "ini_BURN_He_1",
              "BURN_He_75_1",
              "BURN_He_50_1",
              "BURN_He_25_1",
              "end_BURN_He_1",            
              "ini_BURN_C_1",            
              "end_BURN_C_1",
              "ini_BURN_H_2",
              "BURN_H_50_2",
              "BURN_H_25_2",
              "end_BURN_H_2",
              "ini_BURN_He_2",
              "BURN_He_75_2",
              "BURN_He_50_2",
              "BURN_He_25_2",
              "end_BURN_He_2",            
              "ini_BURN_C_2",            
              "end_BURN_C_2",
              "ini_RLOF_A",
              "end_RLOF_A",
              "ini_RLOF_B",
              "end_RLOF_B",
              "ini_RLOF_C",
              "end_RLOF_C",
              "stripping_1",
              "stripping_2",
              "end_1",
              "end_2",
              "unstable_RLOF_A_PA_IV",
              "unstable_RLOF_A_ERK",
              "unstable_RLOF_A_PAULI",
              "unstable_RLOF_A_PABLO",
              "unstable_RLOF_B_PA_IV",
              "unstable_RLOF_B_ERK",
              "unstable_RLOF_B_PAULI",
              "unstable_RLOF_B_PABLO",
              "unstable_RLOF_C_PA_IV",
              "unstable_RLOF_C_ERK",
              "unstable_RLOF_C_PAULI",
              "unstable_RLOF_C_PABLO"]
    ages = Dict()
    data = Dict()
    for elem in what 
        data[elem]=Dict()
    end 

    for stage in stages 
        ages[stage] = f[!, "age_"*which*"_"*stage][model]
        for elem in what 
            data[elem][stage] = f[!, elem*"_"*stage][model]
        end 
    
    end 
    sorted = sort(ages, byvalue=true)
    #print(sorted)
    print("\n")
    print("History for model $model\n")
    print("STAGE                 |  Age (Myr) ")
    for elem in what 
        @printf(" %10s", elem)
    end 
    print("\n")

    for stage in keys(sorted) 
        isnan(sorted[stage]) && continue
        @printf("%22s|", stage)
        @printf("%12.6f", sorted[stage]/1e6)
        for elem in what 
            @printf(" %10.3f", data[elem][stage])
        end 
        @printf("\n")
    end 
    return @sprintf("%4.2f_%4.2f_%4.2f", f[!,"logM"][model],f[!,"logP"][model], f[!,"q"][model])
end



function print_history(f, model, what)
    stages = ["end_BURN_H",
            "ini_BURN_He",
            "end_BURN_He",            
            "ini_BURN_C",            
            "end_BURN_C",
            "stripping",
            "end"]
    ages = Dict()
    data = Dict()
    for elem in what 
        data[elem]=Dict()
    end 

    for stage in stages 
        ages[stage] = f[!, "age_"*stage][model]
        for elem in what 
            data[elem][stage] = f[!, elem*"_"*stage][model]
        end 

    end 
    sorted = sort(ages, byvalue=true)
    #print(sorted)
    print("\n")
    print("History for model $model\n")
    print("STAGE                 |  Age (Myr) ")
    for elem in what 
        @printf(" %10s", elem)
    end 
    print("\n")

    for stage in keys(sorted) 
        isnan(sorted[stage]) && continue
        @printf("%22s|", stage)
        @printf("%12.6f", sorted[stage]/1e6)
        for elem in what 
            @printf(" %10.6f", data[elem][stage])
        end 
        @printf("\n")
    end 

end
