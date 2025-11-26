#!/users/aercolino/.julia/juliaup/julia-1.11.3+0.x64.linux.gnu/bin/julia
include("defs.jl")
using Format

#%% LOAD DIFFERENT Q MODELS  -  NEW CORR FILE! ! ! 


Zs = Dict( 0.010 => "0p5solar",
           0.015 => "0p75solar",
           0.020 => "correctZ",
           0.025 => "1p25solar",
           0.030 => "1p5solar",
           0.035 => "1p75solar",
           0.040 => "twicesolar",
           #0.001 => "smcish"   
)
Ms = range(2.0, 70.0, step=0.5)

master_dir  ="/vol/aibn133/data1/aercolino/GRID_CACHES/drad/single_he_stars/"

for Z in keys(Zs)
    local_dir = master_dir*@sprintf("Z_%.3f/", Z)
    isdir(local_dir) ? nothing : mkdir(local_dir) 

    hal_dir = "/vol/hal/halraid/davidrad/he_stars/single_grids_extraout/"
    hal_subdir = @sprintf("grid_che_helium_single_yoon_++_withN_%s_CC/", Zs[Z])
    

    for M in Ms
        file = local_dir * @sprintf("%.1f.data", M)
        isfile(file) ? continue : nothing 
        
        @printf("READING FILE %.3f/%.1f (%s   %s)\n", Z, M, hal_dir, hal_subdir)
    
        original_file = hal_dir*hal_subdir*@sprintf("/%.1f/LOGS1/history.data", M)
        awk_command = `awk -f /vol/aibn133/data1/aercolino/SOFTWARE/awk_scripts/cleanup_faster.awk $original_file`
        #print( type(awk_command), "  ", awk_command)
        write(file, read(awk_command))
    end
    
    
    
end 






global dummy_checker = []    

for Z in keys(Zs)
    local_dir = master_dir*@sprintf("Z_%.3f/", Z)

    for M in Ms

        fp = nothing 

        @printf("%.3f/%.1f --- LOADING FILE... ", Z, M)
        file = local_dir* @sprintf("%.1f.data", M)
        h = CSV.read(file, header=1, DataFrame,delim=' ',ignorerepeated=true, missingstring="NaN")
        unique!(h, :model_number; keep=:last)


                

        @printf(", SCANNING... ")

        mod  =  h
        models = range(1,length(mod.model_number))
        mass = mod.star_mass
        #center_h1  = mod.center_h1
        #center_he4 = mod.center_he4
        #center_c12 = mod.center_c12
        total_h = mod.total_mass_h1
        R1  = 10 .^ mod.log_R
        dt = 10 .^ mod.log_dt
        end1 = length(mass)
        ZAMS, TAMS, end_Hburn, ini_Heburn, end_Heburn, ini_Cburn, end_Cburn =  evolutionary_checkpoints(h, end1)
            
        @printf(" [HeS] ")
        @printf("He(%5d-%5d) ",  ini_Heburn , end_Heburn)
        @printf("C (%5d-%5d) ",  ini_Cburn, end_Cburn)
            

        CORR_filename = @sprintf("CORR_DRAD_%5.3f.dat", Z)
        @printf("   WRITING TO FILE: %s\n", CORR_filename)
        CORR_file = master_dir*CORR_filename

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

        function hc(column_name)
            return  h[!, column_name]
        end 

        models = range(1,end1,step=1)

        age = hc("star_age")
        dt = 10 .^ hc("log_dt")
        t = maximum(age) .- age
        preCC_1yr = find_nearest(t,1)

        M1 = hc("star_mass")
        M_co_core = hc("c_core_mass")
        M_he_core = hc("he_core_mass")
        M_o_core = hc("o_core_mass")
        M_si_core = hc("si_core_mass")
        M_fe_core = hc("fe_core_mass")
        Mconv = hc("mass_conv_core")

        end_Neburn = find_nearest(M_o_core, 1)
        end_Oburn = find_nearest(M_si_core, 1)
        end_Siburn = find_nearest(M_fe_core, 1)

        xi2p5e = hc("compactness2_5")
        M_tot_h = hc("total_mass_h1")
        M_tot_he = hc("total_mass_he4")
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
        wind_mdot_rate = hc("star_mdot")

        
        Mconv_He =   maximum(Mconv[ini_Heburn : end_Heburn]  ) 
        print(ini_Cburn, " ", end_Cburn, " ")
        Mconv_C  =   ini_Cburn > 0 ? maximum(Mconv[ini_Cburn :  (end_Cburn > 0 ? end_Cburn : end1)]) : NaN 


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
            column["ini_BURN_He"]   =  ini_Heburn 
            column["end_BURN_He"]   =  end_Heburn 
            column["ini_BURN_C"]    =  ini_Cburn  
            column["end_BURN_C"]    =  end_Cburn  
            column["end_BURN_Ne"]   =  end_Neburn  
            column["end_BURN_O"]    =  end_Oburn  
            column["end_BURN_Si"]   =  end_Siburn  
            column["preCC_1yr"]     =  preCC_1yr  
            column["end"]           =  end1 

            for col in keys(column)
                try
                    d[var_name*"_"*col]  = column[col] > 0 ? assigner( variable[column[col]],  typ=typ, quality_check=quality_check ) : assigner( NaN, typ=typ )
                catch E
                    print(E)
                    @printf("\n The error occured with variable %s\n\n", col)
                    throw(ArgumentError("TERMINATING"))

                end
            end

            return d
        end

                    

        datum = OrderedDict()



        datum["M"]  = assigner(M, typ=".3f")
        datum["Z"]   = assigner(Z,     typ=".3f")

        chain_assign_history("age", age, typ = ".12e")
        chain_assign_history("M", M1, typ=".5f", quality_check = true)
        chain_assign_history("M_he_core", M_he_core, typ=".5f")
        chain_assign_history("M_co_core", M_co_core, typ=".5f")
        chain_assign_history("M_o_core", M_o_core, typ=".5f")
        chain_assign_history("M_si_core", M_si_core, typ=".5f")
        chain_assign_history("M_fe_core", M_si_core, typ=".5f")
        chain_assign_history("xi2p5e", xi2p5e, typ=".5f")
        chain_assign_history("M_tot_h", M_tot_h, typ=".5f")
        chain_assign_history("M_tot_he", M_tot_he, typ=".5f")
        chain_assign_history("R", R1, typ=".5f")
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
        datum["Mconv_max_He"] = assigner(Mconv_He,     typ=".5f")
        datum["Mconv_max_C"] = assigner(Mconv_C,     typ=".5f")

        if fp === nothing
            if( ! isfile(CORR_file))
                @printf("CORR file not found. Creating new file.\n")
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
        #@printf("...done\n")

        close(fp)

    end 
        
end

                
                
                
                
                    
                    
                    
