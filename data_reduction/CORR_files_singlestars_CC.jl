#!/users/aercolino/.julia/juliaup/julia-1.10.5+0.x64.linux.gnu/bin/julia
include("defs.jl")
using Format
#%% LOAD DIFFERENT Q MODELS  -  NEW CORR FILE! ! ! 

TO_SAVE = true
master_dir  = "/vol/aibn133/data1/aercolino/GRID_CACHES/hjin/reduced_grid/"

logMs=range(0.700, 2.001, step=0.020)
#logMs=range(1.90, 2.001, step=0.050)
#logPs =range(3.20,3.751, step=0.050)
#logMs=[1.10]
#logPs = [3.40]
#qs_all =  [0.10]
vini = 300

global dummy_checker = []       
fp = nothing 

CORR_filename = @sprintf("CORR_GRID_single_%d.dat", vini)
@printf("--- WRITING TO FILE: %s\n", CORR_filename)
CORR_file = master_dir*CORR_filename


for logM in logMs

        
    local_dir = master_dir*@sprintf("single_stars/")
    
    #READ DATA
    @printf("%.2f  --- LOADING FILES --- ", logM)
    file1 = local_dir*@sprintf("%5.3f_%d.data", logM, vini)
    
    @printf("loading...")
    h = CSV.read(file1, header=1, DataFrame,delim=' ',ignorerepeated=true, missingstring="NaN")
    unique!(h, :model_number; keep=:last)

    @printf(" length: %6d,", length(h.center_h1))

            

    end1  =   length(h.star_age)  
    @printf(" --- SCANNING DATA ---")

    mod  =  h
    center_h1  = mod.center_h1
    center_he4 = mod.center_he4
    center_c12 = mod.center_c12
    total_h = mod.total_mass_h1
    R = 10 .^ mod.log_R         


    ZAMS, TAMS, end_Hburn, ini_Heburn, end_Heburn, ini_Cburn, end_Cburn =  evolutionary_checkpoints(h, end1)
    ix_Rmax  = find_nearest(R, maximum(R))
    stripped  = find_nearest(total_h, 1e-6, out_tol_value = end1, tol = 1e-6  )


    @printf(" [s] ")
    @printf("%7s",  end_Hburn>0 ? " H-dep " : "")
    @printf("%8s",  ini_Heburn>0 ? " He-ign " : "")
    @printf("%8s",  end_Heburn>0 ? " He-dep " : "")
    @printf("%7s",  ini_Cburn>0 ? " C-ign " : "")
    @printf("%7s",  end_Cburn>0 ? " C-dep " : "")
    @printf("  ")
     




    function hc(column_name)
        return h[!, column_name]
    end 
    @printf(  " scanning file... ")
            
    models = range(1,end1,step=1)
    age = hc("star_age")
    dt = 10 .^ hc("log_dt")
    t = maximum(age) .- age
    preCC_1yr = find_nearest(t,1)

    M1 = hc("star_mass")
    M_co_core = hc("c_core_mass")
    M_he_core = hc("he_core_mass")
    M_tot_h = hc("total_mass_h1")
    M_tot_he = hc("total_mass_he4")
    Mconv = hc("mass_conv_core")

    Xc = hc("center_h1")
    Xs = hc("surface_h1")
    Yc = hc("center_he4")
    Ys = hc("surface_he4")
    Cs = hc("surface_c12")
    Cc = hc("center_c12")
    Ns = hc("surface_n14")
    Os = hc("surface_o16")
    Oc = hc("center_o16")
    logL = hc("log_L")
    logTeff = hc("log_Teff")
    log_g = hc("log_g")
    R1  = 10 .^hc("log_R")
    CE_core_R = hc("ce_core_radius")
    Ebind_g = hc("bind_g")
    Ebind_b = hc("bind_b")
    Ebind_h = hc("bind_h")
    Rmax = maxxer(R1)

    lg_wind_mdot_rate = hc("log_abs_mdot")
    v_div_v_crit    =  hc("surf_avg_omega_div_omega_crit")
    v_eq_rot    =  hc("surf_avg_v_rot")
    
    conv_core_top=  min.( hc("conv_mx1_top"), hc("conv_mx2_top")) .* M1
    conv_core_bot=  min.( hc("conv_mx1_bot"), hc("conv_mx2_bot")) .* M1

    Mconv_H  =                   maximum(Mconv[ZAMS       : (      TAMS > 0 ?       TAMS : end1)]) 
    Mconv_He =  ini_Heburn > 0 ? maximum(Mconv[ini_Heburn : (end_Heburn > 0 ? end_Heburn : end1)]) : NaN 
    Mconv_C  =   ini_Cburn > 0 ? maximum(Mconv[ini_Cburn  :       end1]) : NaN 

    function max_core_conv_region(lims, convective_region, topbound, models) 
        conv_top = convective_region[1]
        conv_bot = convective_region[2]
        (lims[1] == 0 || lims[2] == lims[1]) ?  (return NaN) : nothing  
        
        (lims[2] == 0) ? lims[2] = length(topbound) : nothing
        filter_conv_core  = [ (conv_bot[j] < 0.1) && (conv_top[j] < 0.9*topbound[j])   for j in models]
        filter_conv_thick = [ (conv_top[j] .> (conv_bot[j] .+ 0.1) ) for j in models]
        filter_limits     = [ lims[1] <= j <= lims[2] for j in models]
        the_filter        = [ (filter_conv_core[j] && filter_conv_thick[j] && filter_limits[j]) for j in models]
        
        #for j in range(lims[1],lims[2])
        #    @printf("mod %6d | He: %6.3f, Mconv: %.2f, CONV? %5s THICK? %5s IN_LIMITS? %5s\n", 
        #                  j, Yc[j], conv_top[j], filter_conv_core[j], filter_conv_thick[j], filter_limits[j])
        #end 
        #@printf("\nfilter check: %6s => %5.2f\n", true in the_filter, (true in the_filter) ? maximum(conv_top[the_filter]) : NaN)
        
        return (true in the_filter) ? maximum(conv_top[the_filter]) : NaN

    end 


    Mmax_conv_Heburn1 = max_core_conv_region([ini_Heburn, end_Heburn],
                                            [conv_core_top, conv_core_bot], 
                                            M_he_core, 
                                            models) 




    function meanmax(val, lims, dt=dt)
        maximum = -1e99
        (lims[1]==0) && (return [NaN, NaN])

        finish=lims[2]
        (lims[1] > 0 && lims[2] < lims[1]) && (finish = end1)
        limits =  range(lims[1],finish) 
        for j in  limits
            (val[j]>maximum) && (maximum = val[j])
        end
        w_tot = sum(val[limits] .* dt[limits])
        t = sum(dt[limits])
        
        w_mean = t>0 ? (w_tot ./ t) : NaN 
    
        return [w_mean, maximum]
    end



    function chain_assign_history(var_name, variable; d = datum, typ = ".6e", quality_check = true,  secondary=false)
        
        column = OrderedDict()
        column["end_BURN_H"]    =  end_Hburn 
        column["ini_BURN_He"]   = ini_Heburn 
        column["end_BURN_He"]   = end_Heburn 
        column["ini_BURN_C"]    =  ini_Cburn  
        column["end_BURN_C"]    =  end_Cburn  
        column["stripping"]     =   stripped  
        column["end"]           =       end1      

        for col in keys(column)
            try
                d[var_name*"_"*col]  = column[col] > 0 ? assigner( variable[column[col]],  typ=typ, quality_check=quality_check ) : assigner( NaN, typ=typ )
            catch E
                print(E)
                @printf("\n The error occured with variable %s\n\n", col)
                @printf("of value %5d\n", column[col])
                throw(ErrorException("TERMINATING"))

            end
        end

        return d
    end

                
    function chain_assign_mm(var_name, variable; when = ["RLOF", "evo"], d = datum, typ = ".6e")
        if "evo" in when
            res = meanmax(variable, [ZAMS, end_Hburn])
            mean=res[1]; maxim=res[2]
            d[var_name*"_mean_MS"] = assigner(mean, typ=typ)
            d[var_name*"_max_MS"] = assigner(maxim, typ=typ)
            res = meanmax(variable, [ini_Heburn, end_Heburn])
            mean=res[1]; maxim=res[2]
            d[var_name*"_mean_cHeB"] = assigner(mean, typ=typ)
            d[var_name*"_max_cHeB"] = assigner(maxim, typ=typ)
            res = meanmax(variable, [ini_Cburn, end_Cburn])
            mean=res[1]; maxim=res[2]
            d[var_name*"_mean_cCB"] = assigner(mean, typ=typ)
            d[var_name*"_max_cCB"] = assigner(maxim, typ=typ)
        end
        return d
    end

    @printf("storing variables...")

    datum = OrderedDict()



    datum["logM"]  = assigner(logM, typ=".3f")
    
    chain_assign_history("age", age, typ = ".12e")
    chain_assign_history("M", M1, typ=".5f", quality_check = true)
    chain_assign_history("M_he_core", M_he_core, typ=".5f")
    chain_assign_history("M_co_core", M_co_core, typ=".5f")
    chain_assign_history("M_tot_h", M_tot_h, typ=".5f")
    chain_assign_history("M_tot_he", M_tot_he, typ=".5f")
    datum["Mconv_max_Heburn"] = assigner(Mmax_conv_Heburn1,     typ=".5f")
    datum["Mconv_max_H"] = assigner(Mconv_H,     typ=".5f")
    datum["Mconv_max_He"] = assigner(Mconv_He,     typ=".5f")
    datum["Mconv_max_C"] = assigner(Mconv_C,     typ=".5f")

    chain_assign_history("R", R1, typ=".5f")
    chain_assign_history("Rmax", Rmax, typ=".5e")
    chain_assign_history("Xc", Xc, typ=".5f")
    chain_assign_history("Xs", Xs, typ=".5f")
    chain_assign_history("Yc", Yc, typ=".5f")
    chain_assign_history("Ys", Ys, typ=".5f")
    chain_assign_history("Cc", Cc, typ=".5f")
    chain_assign_history("Cs", Cs, typ=".5f")
    chain_assign_history("Ns", Ns, typ=".5f")
    chain_assign_history("Oc", Oc, typ=".5f")
    chain_assign_history("Os", Os, typ=".5f")
    chain_assign_history("logL", logL, typ=".5f")
    chain_assign_history("logT", logTeff, typ=".5f")
    chain_assign_history("log_g", log_g, typ=".5f")
    chain_assign_history("veq_rot", v_eq_rot, typ=".5f")
    chain_assign_history("vcrit_perc", v_div_v_crit, typ=".5f")
    chain_assign_history("CE_core_R", CE_core_R, typ=".5f")
    chain_assign_history("Ebind_g", Ebind_g, typ=".5e")
    chain_assign_history("Ebind_b", Ebind_b, typ=".5e")
    chain_assign_history("Ebind_h", Ebind_h, typ=".5e")


    chain_assign_mm("mdot_wind", 10 .^ lg_wind_mdot_rate)



    if fp === nothing
        if( ! isfile(CORR_file))
            @printf("CORR file not found. Creating new file.")
            fp = open(CORR_file, "w")
            for (key,entry) in datum
                write(fp, cfmt( def_format_width*'s', key ))
            end
            write(fp,"\n")
        else 
            fp = open(CORR_file, "a")
        end
    end

    for (key,entry) in datum
        #print(entry["format"])
        #print(key)
        #print(cfmt(entry["format"], entry["value"]))
        write(fp,  cfmt(entry["format"], entry["value"]))
    end
    write(fp, "\n")   
    @printf("...done\n")

end
@printf("done\n")
close(fp)

            
            
            
            
                
                
                
