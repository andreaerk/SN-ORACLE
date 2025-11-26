#!/users/aercolino/.julia/juliaup/julia-1.11.3+0.x64.linux.gnu/bin/julia
include("defs.jl")
using Format
using HDF5
#%% LOAD DIFFERENT Q MODELS  -  NEW CORR FILE! ! ! 

master_dir  = "/vol/aibn133/data1/aercolino/GRID_CACHES/hjin/reduced_grid/"

logMs=range(0.700, 2.001, step=0.050)
logPs =range(-0.05,3.751, step=0.050)
qs_all =  [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
#logMs=range(0.75, 2.00, step=0.050)
# logPs =range(0.0, .350,step=0.050)
logMs=[2.000]


global dummy_checker = []       

for logM in logMs
    fp = nothing 
    summary_filename = @sprintf("/vol/aibn133/data1/aercolino/GRID_CACHES/hjin/reduced_grid/summary_MW_%5.3f.csv", logM)
    summary_file = CSV.read(summary_filename, header=1, DataFrame)

    for logP in logPs
        
        h = Dict()
        h2 = Dict()
        qs = []
        q_not_finished = []
        terminate_string = Dict()
        # logP = 3.40
        local_dir = master_dir*@sprintf("%5.3f/", logM)
        
        #READ DATA
        @printf("%.2f/____/%.2f  --- LOADING FILES ---\n", logM,logP)
        for q in qs_all
            
            @printf("%.2f/%.2f/%.2f - ", logM,q,logP)
            qt = @sprintf("%.2f", q)
            found = false
            for i in range(1, length(summary_file.q))
                if summary_file.logm[i]==logM && summary_file.logp[i]==logP && summary_file.q[i]==q
                    terminate_string[qt] = summary_file.summary[i]
                    found = true
                    break
                end
            end
            (!found) && (terminate_string[qt] = "Not Run")
            #print(terminate_string[qt])
            try
                terminate_string[qt] = replace(terminate_string[qt], "System reached lower mdot limit<br>" => "")
                terminate_string[qt] = replace(terminate_string[qt], "System reached upper mdot limit<br>" => "")
                terminate_string[qt] = replace(terminate_string[qt], "<br>" => ",")
                terminate_string[qt] = replace(strip(terminate_string[qt]), " " => "_")
            catch
                terminate_string[qt] = "N/A"
            end

            file_id = @sprintf("%5.3f_%5.3f", q, logP)

            is_hdf5  = isfile(local_dir*file_id*"_1_new.hdf5")
            is_ascii = isfile(local_dir*file_id*"_1_new.data")
            h_single = DataFrame() 
            h_sec    = DataFrame() 

            if is_hdf5   
                try
                    @printf("loading hdf5 ...")

                    @printf("loading [1]...")
                    file1 = h5open(local_dir*file_id*"_1_new.hdf5", "r") 
                    for key in keys(file1)
                        h_single[:, key] = read(file1, key)
                    end 

                    @printf("loading [2]...")
                    file2 = h5open(local_dir*file_id*"_2_new.hdf5", "r") 
                    h_sec    = DataFrame() 
                    for key in keys(file2)
                        h_sec[:, key] = read(file2, key)
                    end 

                catch e 
                    print(e)
                    throw(ErrorException)

                    @printf("error reading existing hdf5 files.")
                end
            elseif is_ascii 
                try 
                    file1 = local_dir*file_id*"_1_new.data"
                    file2 = local_dir*file_id*"_2_new.data"
                    
                    @printf("loading [1]...")
                    h_single = CSV.read(file1, header=1, DataFrame,delim=' ',ignorerepeated=true, missingstring="NaN")
                    @printf("loading [2]...")
                    h_sec = CSV.read(file2, header=1, DataFrame,delim=' ',ignorerepeated=true,  missingstring="NaN")

                catch e
                    print(e)
                    throw(ErrorException)
                end 
    
            else 
                h[qt]=[]
                h2[qt]=[]
                @printf("not found. Skipping.\n")
                continue

            end

            try
                dropmissing!(h_single, [:binary_separation])
            catch e 
                h[qt]=[]
                h2[qt]=[]
                @printf("Problem with file entries (missing binary_separation). Skipping.\n")
                continue
            end
            unique!(h_single, :model_number; keep=:last)

            h[qt]=h_single
            unique!(h_sec, :model_number; keep=:last)
            h2[qt] = h_sec
            @printf(" length files: %6d-%6d \n", length(h_single.center_h1), length(h_sec.center_h1))

            append!(qs, q)

            

        end 




        ini_caseA  =   Dict() 
        end_caseA  =   Dict() 
        ini_caseB  =   Dict() 
        end_caseB  =   Dict() 
        ini_caseC  =   Dict() 
        end_caseC  =   Dict() 
        ZAMS       =   Dict() 
        TAMS       =   Dict() 
        H_burn     =   Dict()
        He_burn    =   Dict() 
        H_burn2    =   Dict()
        He_burn2   =   Dict()
        ZAMS2      =   Dict() 
        TAMS2      =   Dict() 
        end_Hburn  =   Dict() 
        ini_Heburn =   Dict() 
        end_Heburn =   Dict() 
        ini_Cburn  =   Dict() 
        end_Cburn  =   Dict() 
        stripped   =   Dict() 
        end_Hburn2 =   Dict() 
        ini_Heburn2=   Dict() 
        end_Heburn2=   Dict() 
        ini_Cburn2 =   Dict() 
        end_Cburn2 =   Dict() 
        stripped2  =   Dict() 
        end1       =   Dict() 
        end2       =   Dict() 
        unstableRLOF_A = Dict() 
        unstableRLOF_B = Dict()
        unstableRLOF_C = Dict()
        STATUS = Dict()

        for q in qs_all
            qt = @sprintf("%.2f", q)
            end1[qt] =  q in qs ?  length(h[qt].star_age) : 0
            end2[qt] =  q in qs ?  length(h2[qt].star_age) : 0
            ! (q in qs) ? continue : nothing
            dummy_checker = [h[qt], h2[qt]]

            ((end1[qt] == end2[qt]) &&  (abs(h[qt].star_age[end1[qt]] - h2[qt].star_age[end2[qt]]) <= 10.) ) ? continue : nothing 

            if abs(end1[qt]-end2[qt]) < 10 ||  (abs(h[qt].star_age[end1[qt]] - h2[qt].star_age[end2[qt]]) <= 10. )
                @printf("%5.2f: CUT (delta age, mod) = %5.3f %2d", q, abs(h[qt].star_age[end1[qt]] - h2[qt].star_age[end2[qt]]),abs(end1[qt]-end2[qt]) )
                new_end = min(end1[qt], end2[qt])
                delete!(h[qt],  [ix for ix in range(new_end+1, end1[qt])])
                delete!(h2[qt], [ix for ix in range(new_end+1, end2[qt])])
                @printf("performed a cut: from [%4d-%4d]",  end1[qt], end2[qt]) 
                end1[qt] = length(h[qt].star_age)
                end2[qt] = length(h2[qt].star_age)
                @printf("to [%4d-%4d]\n", end1[qt], end2[qt]) 

            end
        end 
        @printf("%.2f/____/%.2f  --- SCANNING DATA ---\n", logM,logP)



        for q in qs_all
            qt = @sprintf("%.2f", q)
            @printf("%.2f/%.2f/%.2f ", logM, q, logP)

            if q in  qs

                mod  =  h[qt]
                models = range(1,length(mod.model_number))
                mass = mod.star_mass
                lg_mtransfer_rate = mod.lg_mtransfer_rate
                rl_overflow =  mod.rl_relative_overflow_1
                rl2_overflow = mod.rl_relative_overflow_2
                total_h = mod.total_mass_h1
                upper_mdot = mod.mdot_limit_high
                R1  = 10 .^ mod.log_R
                Rl1 = mod.rl_1
                Rl1out = eval_RLout(mod.star_2_mass ./ mod.star_mass) .* Rl1
                dt = 10 .^ mod.log_dt



                ZAMS[qt], H_burn[qt], TAMS[qt], end_Hburn[qt], ini_Heburn[qt], He_burn[qt], end_Heburn[qt], ini_Cburn[qt], end_Cburn[qt] =  evolutionary_checkpoints(h[qt], end1[qt])
                stripped[qt]  = find_nearest(total_h, 1e-6, out_tol_value = end1[qt], tol = 1e-6  )
                ( stripped[qt]  == end1[qt]) && (  stripped[qt]  = 0)
                
                redo = false
                if ! (end_Hburn[qt] ==0 )
                    for j in range(end_Hburn[qt]+1,  min(end1[qt], end2[qt]) )
        
                        (rl2_overflow[j] <= -0.06) && continue
                        @printf(" IMT! @%5d : %5d-%5d ", j, end1[qt],end2[qt])
                        delete!(h[qt],  [ix for ix in range(j+1, end1[qt])])
                        delete!(h2[qt], [ix for ix in range(j+1, end2[qt])])
                        end1[qt] = length(h[qt].star_age)
                        end2[qt] = length(h2[qt].star_age)
                        @printf(" => %5d-%5d \n", end1[qt],end2[qt])
                        @printf("%15s", "")
                        redo = true
                        break
                    end 
                    if redo 
                        ( TAMS[qt]       > end1[qt]) && (  TAMS[qt] = 0)
                        ( ZAMS[qt]       > end1[qt]) && (  ZAMS[qt] = 0)
                        ( H_burn[qt]["0.50"] > end1[qt]) && ( H_burn[qt]["0.50"] = 0)
                        ( H_burn[qt]["0.25"] > end1[qt]) && ( H_burn[qt]["0.25"] = 0)
                        ( He_burn[qt]["0.75"]> end1[qt]) && ( He_burn[qt]["0.75"] = 0)
                        ( He_burn[qt]["0.50"]> end1[qt]) && ( He_burn[qt]["0.50"] = 0)
                        ( He_burn[qt]["0.25"]> end1[qt]) && ( He_burn[qt]["0.25"] = 0)
                        ( end_Hburn[qt]  > end1[qt]) && (  end_Hburn[qt] = 0)
                        ( ini_Heburn[qt] > end1[qt]) && (  ini_Heburn[qt]= 0)
                        ( end_Heburn[qt] > end1[qt]) && (  end_Heburn[qt]= 0)
                        ( ini_Cburn[qt]  > end1[qt]) && (  ini_Cburn[qt] = 0)
                        ( end_Cburn[qt]  > end1[qt]) && (  end_Cburn[qt] = 0)
                        ( stripped[qt]   > end1[qt]) && (  stripped[qt]  = 0)
                        models = range(1,length(h[qt].model_number))

                    end 
                end

                is_RLOF     = lg_mtransfer_rate .> -7
                Hdep  = ( end_Hburn[qt] == 0) ? end1[qt] :  end_Hburn[qt]
                Hedep = (end_Heburn[qt] == 0) ? end1[qt] : end_Heburn[qt]
                A_ps   =  [ is_RLOF[j] &&         (j <=  Hdep   ) for j in models]
                B_ps   =  [ is_RLOF[j] && ( Hdep < j <= Hedep   ) for j in models]
                C_ps   =  [ is_RLOF[j] && (Hedep < j <= end1[qt]) for j in models]
                ini_caseA[qt] = findfirst(==(true), A_ps)
                (ini_caseA[qt] === nothing) && (ini_caseA[qt] = 0)
                ini_caseB[qt] = findfirst(==(true), B_ps)
                (ini_caseB[qt] === nothing) && (ini_caseB[qt] = 0)
                ini_caseC[qt] = findfirst(==(true), C_ps)
                (ini_caseC[qt] === nothing) && (ini_caseC[qt] = 0)
                end_caseA[qt] = (ini_caseA[qt] > 0) ? findlast(==(true), A_ps) : 0
                (end_caseA[qt] == end1[qt]) && (end_caseA[qt] = 0)
                end_caseB[qt] = (ini_caseB[qt] > 0) ? findlast(==(true), B_ps) : 0
                (end_caseB[qt] == end1[qt]) && (end_caseB[qt] = 0)
                end_caseC[qt] = (ini_caseC[qt] > 0) ? findlast(==(true), C_ps) : 0
                (end_caseC[qt] == end1[qt]) && (end_caseC[qt] = 0)

                function unstable_mass_transfer(lims)
                    #PABLO"s CRITERION: HARD UPPER_MDOT
                    #PAULI"s CRITERION: SOFT UPPER_MDOT
                    #PA_IV"s CRITERION: HARD OUTFLOW
                    #ERK  "s CRITERION: DELAYED OUTFLOW
                    dM_ERK = 0
                    ix_PABLO, ix_PAULI, ix_PA_IV, ix_ERK = 0, 0, 0, 0
                    #Adyn = 0.05
                    
                    finish=lims[2]
                    if lims[1]==lims[2]
                        return Dict("ERK"=>0, "PA_IV"=>0, "PAULI"=>0, "PABLO"=>0)
                    elseif lims[1]>0 && lims[2]==0
                        finish=min(end1[qt], end2[qt])
                    end
                    #orb_p = period[2:end]
                    #dadt = diff(a_sep)./dt[2:end]
                    #dadt_a = dadt./a_sep[2:end]
                    #tau_change = abs.(dadt_a)

    
                    for j in range(lims[1],finish,step=1)
                        if ix_PABLO==0 && lg_mtransfer_rate[j] >= upper_mdot[j]
                            ix_PABLO = j
                        end
                        if R1[j] > Rl1out[j]
                            (ix_PA_IV == 0) ? ix_PA_IV = j : nothing
                            dM_ERK += 10 .^ lg_mtransfer_rate[j] .* dt[j]
                            (ix_ERK == 0 && dM_ERK>= 1) ? ix_ERK = j : nothing
                        end
    
                        #if (tau_change[j-1]*orb_p[j-1]) > Adyn  && !IVA_dyn
                        #    IVA_dyn = true
                        #    ix_IVA = j
                        #end
                    end

                    ix_PAULI =  maximum(lg_mtransfer_rate[lims[1]:finish])>upper_mdot[lims[1]] ? ix_PABLO : 0 

                    return Dict("ERK"=>ix_ERK, "PA_IV"=>ix_PA_IV, "PAULI"=>ix_PAULI, "PABLO"=>ix_PABLO)
                end
                unstableRLOF_A[qt] = unstable_mass_transfer([ini_caseA[qt], end_caseA[qt]])
                unstableRLOF_B[qt] = unstable_mass_transfer([ini_caseB[qt], end_caseB[qt]])
                unstableRLOF_C[qt] = unstable_mass_transfer([ini_caseC[qt], end_caseC[qt]])

                STATUS[qt] = nothing
                (mod.center_he4[end1[qt]] <= 0.01) ? STATUS[qt] = "Ok(HeCC)"  : nothing
                (mod.center_he4[end1[qt]] <= 0.01 && mass[end1[qt]] < 1.38) ? STATUS[qt] = "Ok(HeWD)"  : nothing
                (mod.center_he4[end1[qt]] <  0.01 && mod.center_c12[end1[qt]] <= 0.01) ? STATUS[qt] = "Ok(C-CC)"  : nothing
                (mod.center_he4[end1[qt]] <  0.01 && mod.center_c12[end1[qt]] <= 0.01 && mass[end1[qt]] < 1.38) ? STATUS[qt] = "Ok(COWD)"  : nothing
                (lg_mtransfer_rate[end1[qt]] >= -1) ? STATUS[qt] = "mdot_-1"  : nothing
                (redo) && ( STATUS[qt] = "IMT")
                STATUS[qt] === nothing ? STATUS[qt] = "AbnTerm" : nothing


                mod2 =  h2[qt]
                total_h    = mod2.total_mass_h1
                ZAMS2[qt], H_burn2[qt], TAMS2[qt], end_Hburn2[qt], ini_Heburn2[qt], He_burn2[qt], end_Heburn2[qt], ini_Cburn2[qt], end_Cburn2[qt] =  evolutionary_checkpoints(h2[qt], end2[qt])
                stripped2[qt]  = find_nearest(total_h, 1e-6, out_tol_value = end2[qt], tol = 1e-6  )
                ( stripped2[qt]  == end2[qt]) && (  stripped2[qt]  = 0)



            else
                TAMS[qt]=0
                ZAMS[qt]=0
                H_burn[qt]  =Dict("0.25"=>0, "0.50"=>0)
                He_burn[qt] =Dict("0.25"=>0, "0.50"=>0, "0.75"=>0)
                H_burn2[qt] =Dict("0.25"=>0, "0.50"=>0)
                He_burn2[qt]=Dict("0.25"=>0, "0.50"=>0, "0.75"=>0)
                TAMS2[qt] = 0
                ZAMS2[qt] = 0
                ini_caseA[qt] = 0
                end_caseA[qt] = 0
                ini_caseB[qt] = 0
                end_caseB[qt] = 0
                ini_caseC[qt] = 0
                end_caseC[qt] = 0
                end_Hburn[qt] = 0
                ini_Heburn[qt] = 0
                end_Heburn[qt] = 0
                ini_Cburn[qt] = 0
                end_Cburn[qt] = 0
                stripped[qt] = 0
                end_Hburn2[qt] = 0
                ini_Heburn2[qt] = 0
                end_Heburn2[qt] = 0
                ini_Cburn2[qt] = 0
                end_Cburn2[qt] = 0
                stripped2[qt] = 0
                unstableRLOF_A[qt] = Dict("ERK"=>0, "PA_IV"=>0, "PAULI"=>0, "PABLO"=>0)
                unstableRLOF_B[qt] = Dict("ERK"=>0, "PA_IV"=>0, "PAULI"=>0, "PABLO"=>0)
                unstableRLOF_C[qt] = Dict("ERK"=>0, "PA_IV"=>0, "PAULI"=>0, "PABLO"=>0)
                STATUS[qt] = "n.a." 

            end
    
            #@printf("MT: %6d-%6d  %6d-%6d  %6d-%6d  NUC1: %6d  %6d-%6d  %6d-%6d, str:%6d | %6d | NUC2: %6d  %6d-%6d  %6d-%6d, str:%6d | %6d | \n",
            # ini_caseA[qt], end_caseA[qt], ini_caseB[qt], end_caseB[qt], ini_caseC[qt], end_caseC[qt],
            # end_Hburn[qt], ini_Heburn[qt], end_Heburn[qt], ini_Cburn[qt], end_Cburn[qt], stripped[qt], end1[qt],
            # end_Hburn2[qt], ini_Heburn2[qt], end_Heburn2[qt], ini_Cburn2[qt], end_Cburn2[qt], stripped2[qt], end2[qt])
            
            @printf("  [B] ")
            @printf("A (%5d-%5d) ",  ini_caseA[qt], end_caseA[qt])
            @printf("B (%5d-%5d) ",  ini_caseB[qt] , end_caseB[qt])
            @printf("C (%5d-%5d) ",  ini_caseC[qt], end_caseC[qt])

            
            @printf(" [1] ")
            @printf("H (%5d-%5d) ",  ZAMS[qt], end_Hburn[qt])
            @printf("He(%5d-%5d) ",  ini_Heburn[qt] , end_Heburn[qt])
            @printf("C (%5d-%5d) ",  ini_Cburn[qt], end_Cburn[qt])
            @printf(" [2] ")
            @printf("H (%5d-%5d) ",  ZAMS2[qt], end_Hburn2[qt])
            @printf("He(%5d-%5d) ",  ini_Heburn2[qt] , end_Heburn2[qt])
            @printf("C (%5d-%5d) ",  ini_Cburn2[qt], end_Cburn2[qt])
            @printf("\n")
            
        #SAVE ONTO FILE


        end

     

        CORR_filename = @sprintf("CORR_GRID_%5.3f.dat", logM)
        @printf("%.2f/____/%.2f  --- WRITING TO FILE: %s\n", logM,logP,CORR_filename)
        CORR_file = master_dir*CORR_filename

        function assigner(value; fmt=def_format_width, typ=".6e", quality_check = false)
            ismissing(value) && return Dict("value" => NaN, "format" => fmt*typ)
            isnan(value)     && return Dict("value" => NaN, "format" => fmt*typ)
            
            if quality_check
                if abs(value) >= 1e16 && 'f' in typ
                    value = NaN 
                    typ = ".6e"
                end 
            end 

            return Dict("value" => value, "format" => fmt*typ)
        end 

        for q in qs_all
            #@printf("%.2f/%.2f/%.2f ... ", logM,q, logP)
            qt = @sprintf("%.2f", q)

            function hc(column_name)
                return q in qs ? h[qt][!, column_name] : Vector([NaN])
            end 
            function hc2(column_name)
                return q in qs ? h2[qt][!, column_name] : Vector([NaN])
            end 
            #@printf("%s", q in qs ? " scanning file... " : "setting to nans...")
            
            models = range(1,end1[qt],step=1)
            models2 = range(1,end2[qt],step=1)

            period = hc("period_days")
            age = hc("star_age")
            age2 = hc2("star_age")
            a_sep = hc("binary_separation")
            dt = 10 .^ hc("log_dt")
            t = maximum(age) .- age
            preCC_1yr = find_nearest(t,1)

            M1 = hc("star_mass")
            M2 = hc("star_2_mass")
            Mconv = hc("mass_conv_core")
            Mconv2 = hc2("mass_conv_core")
            M2_2 = hc2("star_mass")
            M_co_core = hc("c_core_mass")
            M_co_core2 = hc2("c_core_mass")
            M_he_core = hc("he_core_mass")
            M_he_core2 = hc2("he_core_mass")
            M_tot_h = hc("total_mass_h1")
            M_tot_h2 = hc2("total_mass_h1")
            M_tot_he = hc("total_mass_he4")
            M_tot_he2 = hc2("total_mass_he4")
            Xc = hc("center_h1")
            Xc2 = hc2("center_h1")
            Xs = hc("surface_h1")
            Xs2 = hc2("surface_h1")
            Yc = hc("center_he4")
            Yc2 = hc2("center_he4")
            Ys = hc("surface_he4")
            Ys2 = hc2("surface_he4")
            Cs = hc("surface_c12")
            Cs2 = hc2("surface_c12")
            Cc = hc("center_c12")
            Cc2 = hc2("center_c12")
            Ns = hc("surface_n14")
            Ns2 = hc2("surface_n14")
            Os = hc("surface_o16")
            Os2 = hc2("surface_o16")
            Oc = hc("center_o16")
            Oc2 = hc2("center_o16")
            logL = hc("log_L")
            logTeff = hc("log_Teff")
            log_g = hc("log_g")
            logL2 = hc2("log_L")
            logTeff2 = hc2("log_Teff")
            log_g2 = hc2("log_g")
            R1  = 10 .^hc("log_R")
            R2 =  10 .^hc2("log_R")
            Rl1 = hc("rl_1")
            Rl2 = hc("rl_2")
            CE_core_R = hc("ce_core_radius")
            CE_core_R2 = hc2("ce_core_radius")
            Ebind_g = hc("bind_g")
            Ebind_b = hc("bind_b")
            Ebind_h = hc("bind_h")
            Ebind_g2 = hc2("bind_g")
            Ebind_b2 = hc2("bind_b")
            Ebind_h2 = hc2("bind_h")
            Rl1out = eval_RLout(M2 ./ M1) .* Rl1
            Rl2out = eval_RLout(M1 ./ M2) .* Rl2

            lg_mtransfer_rate = hc("lg_mtransfer_rate")
            lg_wind_mdot_rate = hc("lg_wind_mdot_1")
            upper_mdot = hc("mdot_limit_high")
            lg_wind_mdot_rate_sec =  hc2("lg_wind_mdot_2")
            v_div_v_crit    =  hc("surf_avg_omega_div_omega_crit")
            v_div_v_crit2    =  hc2("surf_avg_omega_div_omega_crit")
            v_eq_rot    =  hc("surf_avg_v_rot")
            v_eq_rot2    =  hc2("surf_avg_v_rot")
            v_orb_1  = hc("v_orb_1")
            v_orb_2  = hc("v_orb_2")
            
            conv_core_top=  min.( hc("conv_mx1_top"), hc("conv_mx2_top")) .* M1
            conv_core_bot=  min.( hc("conv_mx1_bot"), hc("conv_mx2_bot")) .* M1
            conv_core_top2=  min.( hc2("conv_mx1_top"), hc2("conv_mx2_top")) .* M2_2
            conv_core_bot2=  min.( hc2("conv_mx1_bot"), hc2("conv_mx2_bot")) .* M2_2

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
            #if ini_Heburn[qt]>1
            #     burn_He_ended = end_Heburn[qt]>0 ? end_Heburn[qt] : end1[qt] 
            #     filter_conv_core  =  [ (conv_core_bot[j] < 0.1) && (conv_core_top[j] < 0.9*M_he_core[j])   for j in models]
            #     filter_conv_thick =  [ (conv_core_top[j] .> (conv_core_bot[j] .+ 0.1) ) for j in models]
            #     filter_He_burn    =  [ ini_Heburn[qt] <= j <= burn_He_ended for j in models]
            #     filter_Heburn_conv = [ (filter_conv_core[j] && filter_conv_thick[j] && filter_He_burn[j]) for j in models]
            #     if 1 in filter_Heburn_conv
            #         Mmax_conv_Heburn = maximum(conv_core_top[filter_Heburn_conv])
            #     else
            #         Mmax_conv_Heburn = NaN
            #     end
            # else
            #     Mmax_conv_Heburn = NaN
            # end


            Mmax_conv_Heburn1 = max_core_conv_region([ini_Heburn[qt], end_Heburn[qt]],
                                                    [conv_core_top, conv_core_bot], 
                                                    M_he_core, 
                                                    models) 
            Mmax_conv_Heburn2 = max_core_conv_region([ini_Heburn2[qt], end_Heburn2[qt]],
                                                    [conv_core_top2, conv_core_bot2], 
                                                    M_he_core2, 
                                                    models2) 

            Mconv_H_1 =         ZAMS[qt] > 0 ? maximum(Mconv[ZAMS[qt]       : (       TAMS[qt] > 0 ? TAMS[qt]        : end1[qt]) ]  ) : NaN
            Mconv_He_1 =  ini_Heburn[qt] > 0 ? maximum(Mconv[ini_Heburn[qt] : ( end_Heburn[qt] > 0 ? end_Heburn[qt]  : end1[qt]) ]  ) : NaN
            Mconv_C_1  =   ini_Cburn[qt] > 0 ? maximum(Mconv[ini_Cburn[qt]  :  end1[qt]]) : NaN
            Mconv_H_2 =        ZAMS2[qt] > 0 ? maximum(Mconv2[ZAMS2[qt]      : (      TAMS2[qt] > 0 ? TAMS2[qt]       : end2[qt]) ]  ) : NaN
            Mconv_He_2 = ini_Heburn2[qt] > 0 ? maximum(Mconv2[ini_Heburn2[qt]: (end_Heburn2[qt] > 0 ? end_Heburn2[qt] : end2[qt]) ]  ) : NaN
            Mconv_C_2 =   ini_Cburn2[qt] > 0 ? maximum(Mconv2[ini_Cburn2[qt] : end2[qt]]) : NaN
                                                                                

            function meanmax(val, lims, dt=dt)
                maximum = -1e99
                (lims[1]==0) && (return [NaN, NaN])

                finish=lims[2]
                (lims[1] > 0 && lims[2] < lims[1]) && (finish = end1[qt])
                limits =  range(lims[1],finish) 
                for j in  limits
                    (val[j]>maximum) && (maximum = val[j])
                end
                w_tot = sum(val[limits] .* dt[limits])
                t = sum(dt[limits])
                
                w_mean = t>0 ? (w_tot ./ t) : NaN 
            
                return [w_mean, maximum]
            end

            function eval_deltaM_RLOF(lims)
                deltaM_proper = 0
                finish=lims[2]
                (lims[1]==0) && (return NaN)
                (lims[1] > 0 && lims[2] < lims[1]) && (finish = end1[qt])
                for j in range(lims[1],finish)
                  deltaM_proper += 10^lg_mtransfer_rate[j] .* dt[j]
                end
                return deltaM_proper
            end



            function chain_assign_history(var_name, variable; d = datum, typ = ".6e", quality_check = true,  secondary=false)
                column = OrderedDict()
                column["ini_BURN_H_1"]  =  secondary ? ( ZAMS[qt]            > end2[qt] ? 0 :             ZAMS[qt] )  : ZAMS[qt] 
                column["BURN_H_50_1"]   =  secondary ? ( H_burn[qt]["0.50"]  > end2[qt] ? 0 :   H_burn[qt]["0.50"] )  : H_burn[qt]["0.50"] 
                column["BURN_H_25_1"]   =  secondary ? ( H_burn[qt]["0.25"]  > end2[qt] ? 0 :   H_burn[qt]["0.25"] )  : H_burn[qt]["0.25"]
                column["end_BURN_H_1"]  =  secondary ? ( end_Hburn[qt]       > end2[qt] ? 0 :        end_Hburn[qt] )  : end_Hburn[qt] 
                column["ini_BURN_He_1"] =  secondary ? ( ini_Heburn[qt]      > end2[qt] ? 0 :       ini_Heburn[qt] )  : ini_Heburn[qt] 
                column["BURN_He_75_1"]  =  secondary ? ( He_burn[qt]["0.75"] > end2[qt] ? 0 :  He_burn[qt]["0.75"] )  : He_burn[qt]["0.75"] 
                column["BURN_He_50_1"]  =  secondary ? ( He_burn[qt]["0.50"] > end2[qt] ? 0 :  He_burn[qt]["0.50"] )  : He_burn[qt]["0.50"] 
                column["BURN_He_25_1"]  =  secondary ? ( He_burn[qt]["0.25"] > end2[qt] ? 0 :  He_burn[qt]["0.25"] )  : He_burn[qt]["0.25"]
                column["end_BURN_He_1"] =  secondary ? ( end_Heburn[qt]      > end2[qt] ? 0 :       end_Heburn[qt] )  : end_Heburn[qt] 
                column["ini_BURN_C_1"]  =  secondary ? ( ini_Cburn[qt]       > end2[qt] ? 0 :        ini_Cburn[qt] )  : ini_Cburn[qt]  
                column["end_BURN_C_1"]  =  secondary ? ( end_Cburn[qt]       > end2[qt] ? 0 :        end_Cburn[qt] )  : end_Cburn[qt]  
                column["stripping_1"]   =  secondary ? ( stripped[qt]        > end2[qt] ? 0 :         stripped[qt] )  : stripped[qt]  
                column["end_1"]         =  secondary ? ( end1[qt]            > end2[qt] ? 0 :             end1[qt] )  : end1[qt] 
                column["ini_BURN_H_2"]  = !secondary ? ( ZAMS2[qt]           > end1[qt] ? 0 :            ZAMS2[qt] )  : ZAMS2[qt] 
                column["BURN_H_50_2"]   = !secondary ? ( H_burn2[qt]["0.50"] > end1[qt] ? 0 :  H_burn2[qt]["0.50"] )  : H_burn2[qt]["0.50"] 
                column["BURN_H_25_2"]   = !secondary ? ( H_burn2[qt]["0.25"] > end1[qt] ? 0 :  H_burn2[qt]["0.25"] )  : H_burn2[qt]["0.25"]
                column["end_BURN_H_2"]  = !secondary ? ( end_Hburn2[qt]      > end1[qt] ? 0 :       end_Hburn2[qt] )  : end_Hburn2[qt] 
                column["ini_BURN_He_2"] = !secondary ? ( ini_Heburn2[qt]     > end1[qt] ? 0 :      ini_Heburn2[qt] )  : ini_Heburn2[qt] 
                column["BURN_He_75_2"]  = !secondary ? ( He_burn2[qt]["0.75"]> end1[qt] ? 0 : He_burn2[qt]["0.75"] )  : He_burn2[qt]["0.75"] 
                column["BURN_He_50_2"]  = !secondary ? ( He_burn2[qt]["0.50"]> end1[qt] ? 0 : He_burn2[qt]["0.50"] )  : He_burn2[qt]["0.50"] 
                column["BURN_He_25_2"]  = !secondary ? ( He_burn2[qt]["0.25"]> end1[qt] ? 0 : He_burn2[qt]["0.25"] )  : He_burn2[qt]["0.25"]
                column["end_BURN_He_2"] = !secondary ? ( end_Heburn2[qt]     > end1[qt] ? 0 :      end_Heburn2[qt] )  : end_Heburn2[qt] 
                column["ini_BURN_C_2"]  = !secondary ? ( ini_Cburn2[qt]      > end1[qt] ? 0 :       ini_Cburn2[qt] )  : ini_Cburn2[qt]  
                column["end_BURN_C_2"]  = !secondary ? ( end_Cburn2[qt]      > end1[qt] ? 0 :       end_Cburn2[qt] )  : end_Cburn2[qt]  
                column["stripping_2"]   = !secondary ? ( stripped2[qt]       > end1[qt] ? 0 :        stripped2[qt] )  : stripped2[qt]  
                column["end_2"]         = !secondary ? ( end2[qt]            > end1[qt] ? 0 :             end2[qt] )  : end2[qt] 
                column["ini_RLOF_A"]    = ini_caseA[qt] 
                column["end_RLOF_A"]    = end_caseA[qt] 
                column["end_RLOF_A"]    = end_caseA[qt] 
                column["ini_RLOF_B"]    = ini_caseB[qt] 
                column["end_RLOF_B"]    = end_caseB[qt] 
                column["ini_RLOF_C"]    = ini_caseC[qt] 
                column["end_RLOF_C"]    = end_caseC[qt] 
                for crit in ["PA_IV" , "ERK", "PABLO", "PAULI"]
                    column["unstable_RLOF_A_"*crit]    = unstableRLOF_A[qt][crit]
                end
                for crit in ["PA_IV" , "ERK", "PABLO", "PAULI"]
                    column["unstable_RLOF_B_"*crit]    = unstableRLOF_B[qt][crit]
                end
                for crit in ["PA_IV" , "ERK", "PABLO", "PAULI"]
                    column["unstable_RLOF_C_"*crit]    = unstableRLOF_C[qt][crit]
                end

                for col in keys(column)
                    try
                        d[var_name*"_"*col]  = column[col] > 0 ? assigner( variable[column[col]],  typ=typ, quality_check=quality_check ) : assigner( NaN, typ=typ )
                    catch E
                        print(E)
                        @printf("\n The error occured with variable %s\n\n", var_name*"_"*col)
                        throw(ArgumentError("TERMINATING"))

                    end
                end

                return d
            end

                        
            function chain_assign_mm(var_name, variable; when = ["RLOF", "evo"], d = datum, typ = ".6e")
                if "RLOF" in when
                    res = ini_caseA[qt]>0 ? meanmax(variable, [ini_caseA[qt], end_caseA[qt]]) : [NaN, NaN]
                    mean=res[1]; maxim=res[2]
                    d[var_name*"_mean_RLOF_A"] = assigner(mean, typ=typ)
                    d[var_name*"_max_RLOF_A"] = assigner(maxim, typ=typ)
                    res = ini_caseB[qt]>0 ?  meanmax(variable, [ini_caseB[qt], end_caseB[qt]]) : (NaN, NaN)
                    mean=res[1]; maxim=res[2]
                    d[var_name*"_mean_RLOF_B"] = assigner(mean, typ=typ)
                    d[var_name*"_max_RLOF_B"] = assigner(maxim, typ=typ)
                    res = ini_caseC[qt]>0 ?   meanmax(variable, [ini_caseC[qt], end_caseC[qt]]) : (NaN, NaN)
                    mean=res[1]; maxim=res[2]
                    d[var_name*"_mean_RLOF_C"] = assigner(mean, typ=typ)
                    d[var_name*"_max_RLOF_C"] = assigner(maxim, typ=typ)
                end
                if "evo" in when
                    res = meanmax(variable, [ZAMS[qt], end_Hburn[qt]])
                    mean=res[1]; maxim=res[2]
                    d[var_name*"_mean_MS"] = assigner(mean, typ=typ)
                    d[var_name*"_max_MS"] = assigner(maxim, typ=typ)
                    res = meanmax(variable, [ini_Heburn[qt], end_Heburn[qt]])
                    mean=res[1]; maxim=res[2]
                    d[var_name*"_mean_cHeB"] = assigner(mean, typ=typ)
                    d[var_name*"_max_cHeB"] = assigner(maxim, typ=typ)
                    res = meanmax(variable, [ini_Cburn[qt], end_Cburn[qt]])
                    mean=res[1]; maxim=res[2]
                    d[var_name*"_mean_cCB"] = assigner(mean, typ=typ)
                    d[var_name*"_max_cCB"] = assigner(maxim, typ=typ)
                end
                return d
            end

            deltaM_A_RLOF = eval_deltaM_RLOF([ini_caseA[qt],end_caseA[qt]])
            deltaM_B_RLOF = eval_deltaM_RLOF([ini_caseB[qt],end_caseB[qt]])
            deltaM_C_RLOF = eval_deltaM_RLOF([ini_caseC[qt],end_caseC[qt]])

            #@printf("storing variables...")

            datum = OrderedDict()



            datum["logM"]  = assigner(logM, typ=".3f")
            datum["q"]    = assigner(q,     typ=".3f")
            datum["logP"] = assigner(logP,  typ=".3f")
            datum["STATUS"] = Dict("value" => "\""*STATUS[qt]*"\"", "format" => def_format_width*'s') 

            chain_assign_history("age_1", age, typ = ".12e")
            chain_assign_history("M_1", M1, typ=".5f", quality_check = true)
            chain_assign_history("M_he_core_1", M_he_core, typ=".5f")
            chain_assign_history("M_co_core_1", M_co_core, typ=".5f")
            chain_assign_history("M_tot_h_1", M_tot_h, typ=".5f")
            chain_assign_history("M_tot_he_1", M_tot_he, typ=".5f")
            datum["Mconv_max_Heburn_1"] = assigner(Mmax_conv_Heburn1,     typ=".5f")
            datum["Mconv_max_H_1"] = assigner(Mconv_H_1,     typ=".5f")
            datum["Mconv_max_He_1"] = assigner(Mconv_He_1,     typ=".5f")
            datum["Mconv_max_C_1"] = assigner(Mconv_C_1,     typ=".5f")
            chain_assign_history("R_1", R1, typ=".5f")
            chain_assign_history("Xc_1", Xc, typ=".5f")
            chain_assign_history("Xs_1", Xs, typ=".5f")
            chain_assign_history("Yc_1", Yc, typ=".5f")
            chain_assign_history("Ys_1", Ys, typ=".5f")
            chain_assign_history("Cc_1", Cc, typ=".5f")
            chain_assign_history("Cs_1", Cs, typ=".5f")
            chain_assign_history("Ns_1", Ns, typ=".5f")
            chain_assign_history("Oc_1", Oc, typ=".5f")
            chain_assign_history("Os_1", Os, typ=".5f")
            chain_assign_history("logL_1", logL, typ=".5f")
            chain_assign_history("logT_1", logTeff, typ=".5f")
            chain_assign_history("log_g_1", log_g, typ=".5f")
            chain_assign_history("lg_wind_1", lg_wind_mdot_rate, typ=".5f")            
            chain_assign_history("veq_rot_1", v_eq_rot, typ=".5f")
            chain_assign_history("vcrit_perc_1", v_div_v_crit, typ=".5f")
            chain_assign_history("CE_core_R_1", CE_core_R, typ=".5f")
            chain_assign_history("Ebind_g_1", Ebind_g, typ=".5e")
            chain_assign_history("Ebind_b_1", Ebind_b, typ=".5e")
            chain_assign_history("Ebind_h_1", Ebind_h, typ=".5e")

            chain_assign_history("age_2", age2, typ = ".12e", secondary=true)
            chain_assign_history("M_2", M2_2, typ=".5f", quality_check = true, secondary=true)
            chain_assign_history("M_he_core_2", M_he_core2, typ=".5f", secondary=true)
            chain_assign_history("M_co_core_2", M_co_core2, typ=".5f", secondary=true)
            chain_assign_history("M_tot_h_2", M_tot_h2, typ=".5f", secondary=true)
            chain_assign_history("M_tot_he_2", M_tot_he2, typ=".5f", secondary=true)
            datum["Mconv_max_Heburn_2"] = assigner(Mmax_conv_Heburn2,     typ=".5f")
            datum["Mconv_max_H_2"] = assigner(Mconv_H_2,     typ=".5f")
            datum["Mconv_max_He_2"] = assigner(Mconv_He_2,     typ=".5f")
            datum["Mconv_max_C_2"] = assigner(Mconv_C_2,     typ=".5f")

            chain_assign_history("R_2", R2, typ=".5f", secondary=true)
            chain_assign_history("Xc_2", Xc2, typ=".5f", secondary=true)
            chain_assign_history("Xs_2", Xs2, typ=".5f", secondary=true)
            chain_assign_history("Yc_2", Yc2, typ=".5f", secondary=true)
            chain_assign_history("Ys_2", Ys2, typ=".5f", secondary=true)
            chain_assign_history("Cc_2", Cc2, typ=".5f", secondary=true)
            chain_assign_history("Cs_2", Cs2, typ=".5f", secondary=true)
            chain_assign_history("Ns_2", Ns2, typ=".5f", secondary=true)
            chain_assign_history("Oc_2", Oc2, typ=".5f", secondary=true)
            chain_assign_history("Os_2", Os2, typ=".5f", secondary=true)
            chain_assign_history("logL_2", logL2, typ=".5f", secondary=true)
            chain_assign_history("logT_2", logTeff2, typ=".5f", secondary=true)
            chain_assign_history("log_g_2", log_g2, typ=".5f", secondary=true)
            chain_assign_history("lg_wind_2", lg_wind_mdot_rate_sec, typ=".5f", secondary=true)            
            chain_assign_history("veq_rot_2", v_eq_rot2, typ=".5f", secondary=true)
            chain_assign_history("vcrit_perc_2", v_div_v_crit2, typ=".5f", secondary=true)
            chain_assign_history("CE_core_R_2", CE_core_R2, typ=".5f", secondary=true)
            chain_assign_history("Ebind_g_2", Ebind_g2, typ=".5e", secondary=true)
            chain_assign_history("Ebind_b_2", Ebind_b2, typ=".5e", secondary=true)
            chain_assign_history("Ebind_h_2", Ebind_h2, typ=".5e", secondary=true)


            chain_assign_history("a", a_sep, typ=".5f", quality_check = true)
            chain_assign_history("period", period, typ=".5f", quality_check = true)
            chain_assign_history("Rl_1", Rl1, typ=".5f", quality_check = true)
            chain_assign_history("Rl_2", Rl2, typ=".5f", quality_check = true)
            chain_assign_history("Rlout_1", Rl1out, typ=".5f")
            chain_assign_history("Rlout_2", Rl2out, typ=".5f")
            chain_assign_history("v_orb_1", v_orb_1, typ=".5f")
            chain_assign_history("v_orb_2", v_orb_2, typ=".5f")





            
            chain_assign_mm("mdot_RLOF", 10 .^ lg_mtransfer_rate, when = ["RLOF"])
            chain_assign_mm("mdot_wind", 10 .^ lg_wind_mdot_rate)
            datum["deltaM_A"] = assigner(deltaM_A_RLOF, typ = ".5f")
            datum["deltaM_B"] = assigner(deltaM_B_RLOF, typ = ".5f")
            datum["deltaM_C"] = assigner(deltaM_C_RLOF, typ = ".5f")

            datum["summary"] = Dict("value" => "\""*terminate_string[qt]*"\"", "format" => "%200s") 

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
            #@printf("...done\n")

        end
        @printf("%.2f/____/%.2f  File Written \n\n", logM,logP)

    end
    @printf("done\n")
    close(fp)

end 
    
@printf("Building Master-Corr file\n")
master_corr = CSV.read(master_dir*"CORR_GRID_0.700.dat", header=1, DataFrame,delim=' ',ignorerepeated=true)
for logM in 0.750:0.050:2.000 
    _fn = @sprintf("CORR_GRID_%5.3f.dat", logM)
    @printf("inclduing %s\n", _fn)
    _f = CSV.read(master_dir*_fn, header=1, DataFrame,delim=' ',ignorerepeated=true)
    master_corr = vcat(master_corr, _f)
end 

@printf("reducing to hdf5...")
file_hdf5 = h5open(master_dir * "master_CORR.hdf5", "w") 
for col in names(master_corr) 
    file_hdf5[col] = master_corr[!, col]
end   
close(file_hdf5) 
@printf(" done!\n") 

            
            
            
                
                
                
