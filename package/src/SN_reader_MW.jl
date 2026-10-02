# Copyright (C) 2026  [Andrea Ercolino]
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.




function read_hdf5_into_dataframe(filename)
    @printf("Reading Master CORR file %s ...  ", filename)
    a = h5open(filename, "r") do file 
        # Each key in the root is a column
        colnames = keys(file)
    
        # Read each dataset into a NamedTuple
        data = (; (Symbol(name) => read(file[name]) for name in colnames)...)
    
        @printf("done\n")

        return DataFrame(data)
    end
end 

function modkeys(; logM::Union{Float64,String}=-Inf,
                   logP::Union{Float64,String}=-Inf,
                   q::Union{Float64,String}=-Inf)
    if logM isa Float64
        return @sprintf("%4.2f_%4.2f_%4.2f", logM, logP, q)
    else
        return logM * "_" * logP * "_" * q
    end
end

function get_from_modkey(;what::String, key::String=key, returntype=Float64)
    if what=="logM"
        ixs=1:4
    elseif what == "logP"
        ixs=6:9
    elseif what == "q"
        ixs=11:14
    else
        throw(ErrorException("(get_from_modkey) 'what' is not recognized"))
    end

    if returntype==Float64
        return parse(returntype, key[ixs])
    elseif returntype==String
        return key[ixs]
    else
        throw(ErrorException("(get_from_modkey) 'returntype' is not recognized"))
    end
end

function predict_SN_from_each_model(;
        f = nothing,    
        MERGER_CRITERION=nothing,  # X, PA_IV, ERK, PABLO, PAULI
        MERGER_EOL=nothing,        # "TOTAL"   "CORE"   "ANALYTIC"
        MERGER_MASSLOSS=nothing,   # "FIX"+fraction   "Energy"
        EXP_CRIT=nothing,          # X, MM+cal, Ertl+cal, comp_value
        KICKS=nothing,             # 0, Inf, or Distr
        logMs = range(0.70, 2.001, step=0.050),
        logPs = range(0.05, 3.750, step=0.050),
        qs    = range(0.10, 0.951, step=0.050)
        )
    global wie_warning = 0
    global counter_call_gefibH = 0

    f_def = occursin("FIX", MERGER_MASSLOSS) ? parse(Float64, MERGER_MASSLOSS[4:end]) : 0  


    MODELS = Dict()
    for logM in logMs 
        for logP in logPs 
            for q in qs
                MODELS[modkeys(logM=logM,logP=logP,q=q)] = Dict()
            end
        end
    end

    outcomes = [:WD, :IIP, :SN87A, :IIb, :Ibc, :BH, :IIn, :Ibn, :X]
    logM = NaN
    @printf("Scanning logMs\n")

    for model in range(1, length(f.q))
        # (f.logM[model] != logM) && @printf("%.2f  \r", f.logM[model])
        ((model%500)==0) && @printf("%.2f [%5d/%5d] \r", f.logM[model], model, length(f.q)) 
        logM = f.logM[model]
        logP = f.logP[model]
        q = f.q[model]
        case = ""
        terminated1 = ""
        terminated2 = ""
        single = false 
        merger = nothing 

        if ! ( (logM in logMs) && (logP in logPs) && (q in qs))
            continue 
        end 
        
        entry(what, when, which) = f[!, what*"_"*which*"_"*when*"_"*which][model]
        entry(what, when)        = f[!, what*"_"*when][model]

        casus = nothing

        postkick_states = [(0.,0.,0.)]

        if not_run(f, model)
            if f.logP[model] > 3.000
                terminated1 = "not run: single"
                terminated2 = "not run: single"
                single = true
                case = "X"
            elseif f.logP[model] < 3.150 #!!!#
                terminated1 = "not run: ZAMS merger"
                terminated2 = "not run: ZAMS merger"
                case = "A"
            else
                @printf("MODEL NOT RUN BUT AT AN INTERMEDIATE PERIOD!!!\n")
                @printf("logM/logP/q = %.2f/%.2f/%.2f", logM, logP, q)
                print_history(f, model, ["M_1", "M_2"], "2")
                throw(ErrorException("model not run but at an intermediate period? Check!")) 
            end 
        else 
            terminated1 = termination_stage(f, model, "1")
            terminated2 = termination_stage(f, model, "2")
        end 

        merger = is_merger(f, model, MERGER_CRITERION, terminated1, terminated2)
        case = merger["RLOFcase"]

        print_mod_det(bool) = bool ? @printf("model %5d: %5.2f/%5.2f/%.2f = " ,model, logM, logP, q) : nothing
        SN1s = []
        SN2s = []
        MM = Dict()
        t_SN1 = Float32[]
        t_SN2 = Float32[] 
        first_SN = nothing 
        solverflag1=nothing 
        solverflag2=nothing
        v_kicks = Float32[]
        extraMTs = [] 
        stats_ECI = [] 
        post_SN_orbits = []
        pre_SN_orbit = nothing

        SN = nothing
        unboundstate=Dict("ix"=>0, "counter"=>0, "v_kick"=>[])

        if single
            M1 =  10^logM 
            M2 =  10^logM * q
            casus = "Effectively Single"
            extraMTs = [""] 
            v_kicks = Float32[0]
            terminated1 = "not run"
            terminated2 = "not run"
            first_SN = "1"
            second_SN  = "2"
            t_SN1 = [f_age_CB(M1)]
            t_SN2 = [f_age_CB(M2)]
            MM["facc_A"]  = NaN
            MM["facc_B"]  = NaN

            SN, endvals=SINGLE_SN(M1,  EXP_CRIT, "single")
            SN1s=[SN]
            MM["M1max"]    = Float16[M1]
            MM["endvals1"] = [endvals]
            solverflag1 = "sin"
            solverflag2 = "sin"

            eci = ECI.estimate_ECI(;  Mej  = endvals.M_end - endvals.M_remnant_b, #Msun
                                        Eexp = endvals.E_exp, #erg
                                        a  = +Inf, #Rsun
                                        M2 = M2, #Msun 
                                        R2 = 0., #Rsun
                                        RL2= +Inf#Rsun
                                        )
            post_SN_orbit = (a=-Inf, a_peri=-Inf, P=-Inf, e=-Inf, ϵ=-Inf, 
                            orbit = :unknown, 
                            will_one_periastron_occur=false, 
                            when_will_periastron_occur = NaN)
            push!(post_SN_orbits, post_SN_orbit)
            push!(stats_ECI, eci)

            SN2, endvals=SINGLE_SN(M2,  EXP_CRIT, "single")
            MM["M2max"]    =  Float16[M2]
            MM["dM_C"]     =  Float16[0]
            MM["dM_C2"]     = Float16[0]
            MM["endvals2"] = [endvals]
            SN2s = [SN2]
            pre_SN_orbit = (a=+Inf, v1=NaN, v2=NaN)
        elseif merger["merger?"]
            M1 = merger["M1"]
            M2 = merger["M2"]
            M1core = merger["M1core"]
            M2core = merger["M2core"]
            M_core = M1core+M2core
            M1h = merger["M1h"]
            M2h = merger["M2h"]
            Xc = merger["Xc"]
            factor = unstable_mass_transfer_analysis(f,model,MERGER_CRITERION, MERGER_MASSLOSS; val = f_def)
            factor = min(factor,1)
            M = (1-factor)*(M1 + M2)
            first_SN = "M"
            f_acc = (M2)/(10^logM)+(M1-10^logM)/(10^logM)
            second_SN  = "X"
            endvals = Dict()
            extraMTs = [""] 
            v_kicks = Float32[0]
            M_max_pMS_single = interpol_single_Mcore_pMS["M_pMS"](M_core)
            M_max_single = interpol_single_Mcore_pMS["M_Max"](M_core)

            solverflag2 = "X"

            eci = ECI.estimate_ECI(;  Mej  = NaN, #Msun
                                        Eexp = NaN, #erg
                                        a  = NaN, #Rsun
                                        M2 = NaN, #Msun 
                                        R2 = NaN, #Rsun
                                        RL2= NaN#Rsun
                                        )
            push!(stats_ECI, eci)

            if MERGER_EOL == "TOTAL" || occursin("A",  merger["case"]) || (M_core < 0.1 && M2core<0.1)
                solverflag1 = "(rejuvinated) SINGLE_SN"

                SN,endvals["1"]=SINGLE_SN(M,  EXP_CRIT,
                                        occursin("A",  merger["case"]) ? "single" : "accretor", 
                                        f_acc = occursin("A",  merger["case"]) ? NaN : f_acc,
                                        caseC=occursin("IMT-C",  merger["case"]), solver_string=solverflag1)
                
            elseif MERGER_EOL == "CORE" && !(M2core > 0.1 && M1core > 0.1) && !occursin("IMT-C",  merger["case"])
                solverflag1 = "(CORE) SINGLE_SN"

                SN, endvals["1"]=SINGLE_SN(M,  EXP_CRIT, "accretor", 
                                    f_acc = f_acc, 
                                    Mi=10^logM, 
                                    early_or_late=merger["early_or_late"],
                                    caseC=occursin("IMT-C",  merger["case"]), solver_string=solverflag1)

            elseif M2core > 0.1 && M1core > 0.1 && !occursin("IMT-C",  merger["case"])

                f_acc = (M-M_max_single)/(M_max_single)

                if f_acc > 1
                    solverflag1 = "(IMT) SINGLE_SN"
                    SN,endvals["1"]=SINGLE_SN(M,  EXP_CRIT, "accretor", 
                                f_acc = f_acc, 
                                Mi=M_max_single, 
                                early_or_late="late",
                                caseC=occursin("IMT-C",  merger["case"]), solver_string=solverflag1)
                    
                else 
                    solverflag1 = "(IMT) get_end_from_ini_BURN_He"
                    SN, endvals["1"] = get_end_from_ini_BURN_He(M_core, M-M_core, EXP_CRIT; reference_grid="MW-ZAMS", caseC=occursin("IMT-C",  merger["case"]), solver_string=solverflag1)

                end
            elseif occursin("IMT-C",  merger["case"])
                solverflag1 = "(IMT-C) SN_OUTPUT_manual"

                SN, endvals["1"] =SN_OUTPUT_manual(M, M_core, f.M_co_core_1_end_1[model]+f.M_co_core_2_end_1[model], M_max_pMS_single, M1h+M2h, Xc, "end_BURN_He", EXP_CRIT;
                                                deltaM_C = Inf, reference_grid = "MW-ZAMS", solver_string=solverflag1)
            else
                throw(ErrorException)
            end


            SN2 = :X 
            endvals["2"]=output_run_data( )
            SN1s=[SN]

            t_end = f[!, "age_1_end_1"][model]
            SN2s = [:X]
            MM["endvals1"]  = [endvals["1"]]
            MM["endvals2"]  = [endvals["2"]]
            MM["M1max"]     = Float16[M]
            MM["M2max"]     = Float16[NaN]
            MM["facc_A"]    = occursin("A",  merger["case"]) ? f_acc : NaN
            MM["facc_B"]    = occursin("B",  merger["case"]) ? f_acc : NaN
            MM["dM_C"]      = Float16[NaN]
            MM["dM_C2"]     = Float16[NaN]
            pre_SN_orbit = (a=NaN, v1=NaN, v2=NaN)

            post_SN_orbit = (a=NaN, a_peri=NaN, P=NaN, e=NaN, ϵ=NaN, 
                orbit = :merged, 
                will_one_periastron_occur=false, 
                when_will_periastron_occur = NaN)
            push!(post_SN_orbits, post_SN_orbit)

            t_SN1 = Float32[NaN] 
            t_SN2 = Float32[NaN]

            (merger["case"] === nothing) && (print(merger))
            if occursin("B", merger["case"])
                if MERGER_EOL == "TOTAL"
                    # t_SN1 = t_preCC_pMS_M(M)+t_end
                elseif MERGER_EOL == "CORE"
                    # t_SN1 = t_preCC_pMS_Mcore(Mcore)+t_end
                    # t_SN1 = NaN
                end
            elseif occursin("A", merger["case"])
                # t_SN1 = retrieve_tpreCC_mergers([M], [Xc])+t_end
            elseif occursin("C", merger["case"])
                # t_SN1 = t_preCC_pHeB_Mcore(M_core)+t_end
            elseif occursin("IMT", merger["case"])
                if terminated1 in ["pMS", "HeB", "pHeB", "CB", "pCB", "Hedep"] || terminated2 in ["pMS", "HeB", "pHeB", "CB", "pCB", "Hedep"]
                    # t_SN1 = t_preCC_pHeB_Mcore(M_core)+t_end
                    case *= "IMT"
                else
                    print(merger)
                    print(terminated1, terminated2)
                    throw(ErrorException)
                end
            end
            
            casus = merger["case"] * "-" * merger["early_or_late"]
        else 
            M1 = f.M_1_end_1[model]
            M2 = f.M_2_end_2[model]

            #Mcore_he_end =  f.M_he_core_end[model]
            #Mcore_co_end =  f.M_co_core_end[model]
            #Mh = f.M_tot_h_end[model]
            age1 = entry("age_1", "end_1")
            age2 = entry("age_2", "end_2")
            first_SN = (age1 <= age2+10) ? "1" : "2"
            #(first_SN == "2") ? @printf("first SN: %s; deltaT = %.3e\n", first_SN, age2-age1) : nothing 
            second_SN = (first_SN == "1") ? "2" : "1"
            endvals=Dict()
            solverflag1 = "bin1-"
            SN, endvals["1"]=SN_OUTPUT_firstSN(first_SN=="1" ? terminated1 : terminated2, f, model, EXP_CRIT, solverflag1; which=first_SN)
            
            SN1s=[SN]
            t_SN1 = [(first_SN == "1" && (SN!="WD" && SN!="BH")) ? age1 : age2]
            MM["dM_C"]     = Float16[endvals["1"].deltaM_C]
            MM["endvals1"] = [endvals["1"]]
            MM["M1max"]    = first_SN == "1" ? Float16[10 ^ f[!, "logM"][model]] : Float16[ 10 ^ f[!, "logM"][model] *   f[!, "q"][model]]
            solverflag2 = "bin2-"
            MM["facc_A"]  = (f[!, "M_2_end_RLOF_A"][model]-f[!, "M_2_ini_RLOF_A"][model])/ 10^logM*q
            MM["facc_B"]  = (f[!, "M_2_end_RLOF_B"][model]-f[!, "M_2_ini_RLOF_B"][model])/ 10^logM*q

            MM["M2max"]  = Float16[]
            MM["dM_C2"]  = Float16[]
            MM["endvals2"]           = []
            stats_ECI                = []
            SN2s                     = []


            M2_atCC   = entry("M_"*second_SN,    "end_"*first_SN)
            M2_hedep  = entry("M_"*second_SN,    "end_BURN_He_"*first_SN)
            M1_atCC   = entry("M_"*first_SN,     "end_"*first_SN)
            M1_hedep  = entry("M_"*first_SN,     "end_BURN_He_"*first_SN)
            RL2_atCC  = entry("Rl_"*second_SN,   "end_"*first_SN)
            RL2_hedep = entry("Rl_"*second_SN,   "end_BURN_He_"*first_SN)
            R2_atCC  = entry("R_"*second_SN,     "end_"*first_SN)
            R2_hedep = entry("R_"*second_SN,     "end_BURN_He_"*first_SN)
            v1_atCC  = entry("v_orb_"*first_SN,  "end_"*first_SN)
            v1_hedep = entry("v_orb_"*first_SN,  "end_BURN_He_"*first_SN)
            v2_atCC  = entry("v_orb_"*second_SN, "end_"*first_SN)
            v2_hedep = entry("v_orb_"*second_SN, "end_BURN_He_"*first_SN)
            a_atCC   = entry("a",                "end_"*first_SN)
            a_hedep  = entry("a",                "end_BURN_He_"*first_SN)

            pre_SN_orbit =  if SN == :WD 
                                derive_a_massloss_smooth(a_atCC,  M1_atCC,  M2_atCC,  M1_atCC-endvals["1"].M_end)
                            elseif SN == :Ibn 
                                derive_a_massloss_smooth(a_hedep, M1_hedep, M2_hedep, M1_hedep-endvals["1"].M_end)
                            else
                                (a = a_atCC, v1 = v1_atCC, v2 = v2_atCC)
                            end

            kick_case = ! (occursin("A", case) || occursin("B", case)) ? "ZAMS" : ((endvals["1"].deltaM_C > 0 && SN == :Ibn) ? "CaseBB" : "stripped")
            kick_parameters = get_kick_parameter(KICKS, endvals["1"], kick_case)
            (kick_case != "CaseBB" && SN == :Ibn) && (dm = log10(endvals["1"].deltaM_C); st = endvals["1"].solver_string; throw(ErrorException("Ibn not classified as Case BB system? log deltaMC = $dm. Where ? $st")))
            postkick_states = SN in [:WD, :BH] ? [(0.,0.,0.)] : get_kicks(KICKS, kick_parameters, SN in [:WD, :BH] )
            for k in range(1,length(postkick_states))
                kick = postkick_states[k]
                
                post_SN_orbit=nothing
                if ! (SN in [:WD, :BH])
                    v_1_new, v_2_new = apply_kick_and_recenter(endvals["1"].M_remnant_g, pre_SN_orbit.v1, M2_atCC, pre_SN_orbit.v2; v_kick=kick)
                    post_SN_orbit   =  derive_orbital_solution(endvals["1"].M_remnant_g*Msun, 
                                                         [0.,0.,0.]*Rsun, 
                                                         v_1_new*km, 
                                                         M2_atCC*Msun, 
                                                         [pre_SN_orbit.a,0.,0.]*Rsun, 
                                                         v_2_new*km,
                                                         R2_atCC*Rsun+15*km)
                else
                    post_SN_orbit = (a=pre_SN_orbit.a, a_peri=pre_SN_orbit.a, P=NaN, e=NaN, ϵ=NaN, 
                                    orbit = :bound, 
                                    will_one_periastron_occur=false, 
                                    when_will_periastron_occur = NaN)
                end

                eci = nothing
                new_RL = nothing
                if post_SN_orbit.will_one_periastron_occur
                    new_RL = eval_RL(f[!, "M_"*second_SN*"_end_"*first_SN][model], endvals["1"].M_remnant_g, "a", post_SN_orbit.a_peri/Rsun)
                else 
                    new_RL = eval_RL(f[!, "M_"*second_SN*"_end_"*first_SN][model], endvals["1"].M_remnant_g, "a", pre_SN_orbit.a)
                    if post_SN_orbit.when_will_periastron_occur < 0 
                        post_SN_orbit = (a=post_SN_orbit.a, a_peri=post_SN_orbit.a_peri, P=post_SN_orbit.P, e=post_SN_orbit.e, ϵ=post_SN_orbit.ϵ, 
                                       orbit = post_SN_orbit.orbit, 
                                       will_one_periastron_occur=false, 
                                       when_will_periastron_occur = 0.)
                    end
                end
                
                
                eci = ECI.estimate_ECI(;  Mej  = endvals["1"].M_end - endvals["1"].M_remnant_b, #Msun
                                        Eexp = endvals["1"].E_exp, #erg
                                        a    = pre_SN_orbit.a, #Rsun
                                        M2   = M2_atCC, #Msun 
                                        R2   = R2_atCC, #Rsun
                                        RL2  = post_SN_orbit.a_peri/Rsun #new_RL #Rsun
                                        )

                if !eci["ogata"]["ECI?"] && post_SN_orbit.orbit == :unbound 
                    unboundstate["counter"]+=1
                    push!(unboundstate["v_kick"], kick)
                    if unboundstate["ix"] == 0
                        unboundstate["ix"] = k
                    else 
                        continue
                    end
                end

                if (post_SN_orbit.orbit == :unbound &&  SN in [:WD, :BH]) 
                    @printf("WARNING! UNBOUND ORBIT FOLLOWING BH/WD FORMATION! IMPOSSIBLE!\n")
                    @printf("logM-logP-q : %5.3f/%5.3f/%5.3f\n", logM, logP, q)
                    @printf("SN : %s - endvals = ", SN, )
                    print(endvals["1"], "\n")
                    @printf("v_old: %.2f/%.2f; v_new: %.2f/%.2f\n", v1_atCC, v2_atCC, norm(v_1_new), norm(v_2_new))
                    print("kick was ", kick, "\n")
                    throw(ErrorException)
                    
                end
                
                if post_SN_orbit.orbit == :direct_collision
                    SN2 = :X 
                    endvals["2"] = output_run_data( )
                    extraMT = ""
                    casus = "end1+directcollision"
                    t_sn2 = NaN
                else
                    SN2, endvals["2"], extraMT =SN_OUTPUT_secondSN(first_SN=="1" ? terminated2 : terminated2, f, model, EXP_CRIT, norm(kick), second_SN, SN, endvals["1"], pre_SN_orbit, post_SN_orbit, solverflag2)
                    # (SN2["Ibc-i"]>0) && print_history(f, model, ["M_1", "M_he_core_1", "M_2", "M_he_core_2", "R_2", "Rl_2"], "2") 
                    t_sn2 = (second_SN == "1" && (SN2!="WD" && SN2!="BH")) ? age1 : age2
                    casus = "end"
                end
                push!(v_kicks, norm(kick))
                push!(post_SN_orbits, post_SN_orbit)
                push!(t_SN2, t_sn2)
                push!(extraMTs, extraMT)
                push!(SN2s, SN2)
                push!(stats_ECI, eci)
                push!(MM["dM_C2"], NaN)
                push!(MM["M2max"], first_SN == "2" ?  10 ^ f[!, "logM"][model] *   f[!, "q"][model] : 10 ^ f[!, "logM"][model])
                push!(MM["endvals2"], endvals["2"]        )
            end

            if (length(postkick_states) != (length(v_kicks)+max(unboundstate["counter"],1)-1))
                @printf("Warning!\n Postkick states generated = %5d\n", length(postkick_states))
                @printf("kick states returned %5d, of which %5d are unbound states  ", length(v_kicks), unboundstate["counter"])
                throw(ErrorException("generated !+ returned + (unbound - 1)"))
            end
            
        end

        (nothing in SN1s) &&  (throw(ErrorException("SN 1 has nothing")))
        (nothing in SN2s) &&  (throw(ErrorException("SN 2 has nothing")))

        modkey=modkeys(logM=logM,logP=logP,q=q)
        MODELS[modkey]["mod"]    = model
        MODELS[modkey]["1"] = (SN      = copy(SN1s),
                                            Mmax    = MM["M1max"],
                                            dM_C    = MM["dM_C"],
                                            endvals = MM["endvals1"],
                                            t_SN    = t_SN1,
                                            case    = [case]
                                            )
        MODELS[modkey]["2"] = (SN      = copy(SN2s),
                                            Mmax    = MM["M2max"],
                                            dM_C    = MM["dM_C2"],
                                            endvals = MM["endvals2"],
                                            t_SN    = t_SN2,
                                            case    = extraMTs
                                            )
        MODELS[modkey]["1_whichstar"] = first_SN
        # MODELS[modkey]["facc_A"] = MM["facc_A"] 
        # MODELS[modkey]["facc_B"] = MM["facc_B"] 
        MODELS[modkey]["1st_kick_whole"] = postkick_states
        MODELS[modkey]["unboundstate"] = unboundstate
        MODELS[modkey]["1st_kick"] = v_kicks
        MODELS[modkey]["casus"] = casus
        MODELS[modkey]["end_primary"]   = terminated1
        MODELS[modkey]["end_secondary"] = terminated2
        MODELS[modkey]["ECI"] = stats_ECI
        MODELS[modkey]["post_SN_orbit"] = post_SN_orbits
        MODELS[modkey]["pre_SN_orbit"] = pre_SN_orbit
        
    end
    @printf("\n")
    wie_warning > 0 && @printf("Warning: will_it_explode had to derive properties outside the safe monotonic regions %4d times\n", wie_warning)
    counter_call_gefibH > 0 && @printf("Warning: get_end_from_ini_BURN_He called %4d times\n", counter_call_gefibH)
    return MODELS 

end



SANA_logP_pdf(x, power) = max(x, 0.05) ^ (power)
SANA_q_pdf(x,    power) = x ^ (power)
Salpeter_IMF(x,  power) = x ^ (power)
function integrate_logM(logMmin, logMmax, power)
    m_min = 10 .^ logMmin
    m_max = 10 .^ logMmax
    return quadgk(m -> Salpeter_IMF(m, power),  m_min, m_max, rtol=1e-9)
end 
function integrate_logM_convolve(logMmin, logMmax, power, f = x -> x* 0 )
    m_min = 10 .^ logMmin
    m_max = 10 .^ logMmax
    return quadgk(m -> Salpeter_IMF(m, power)*f(m),  m_min, m_max, rtol=1e-9)
end 



#ASSIGN A INITIAL-PROBABILITY DISTRIBUTION TO EACH SYSTEM 
function do_SN_popsynth(;
    f = nothing,
    MERGER_CRITERION = nothing, 
    MERGER_EOL = nothing, 
    MERGER_MASSLOSS = nothing, 
    EXP_CRIT = nothing, 
    KICKS = nothing,
    pow_m = nothing, 
    pow_p = nothing, 
    pow_q = nothing, 
    fB = nothing, 
    fB_value = nothing,
    logMs = range(0.70, 2.001, step=0.050),
    logPs = range(0.05, 3.70, step=0.05),
    qs    = range(0.10, 0.951, step=0.050),
    WRITE_TO_FILE = true,
    MODELS = nothing,
    only_SN = false
    )


    if MODELS === nothing
        MODELS = predict_SN_from_each_model(; 
                f = f,
                MERGER_CRITERION = MERGER_CRITERION,
                MERGER_EOL = MERGER_EOL,
                MERGER_MASSLOSS = MERGER_MASSLOSS, 
                EXP_CRIT = EXP_CRIT, 
                KICKS = KICKS,
                logMs = logMs,
                logPs = logPs,
                qs    = qs,
                )
    end 

    outcomes = [:WD, :IIP, :SN87A, :IIb, :Ibc, :IIn, :Ibn, :BH, :X]
    SN_outcomes = copy(outcomes)
    filter!(!=(:WD),SN_outcomes)
    filter!(!=(:X), SN_outcomes)
    filter!(!=(:BH),SN_outcomes)

    delta = 0.050/2

    MERGER_MASSLOSS *=  MERGER_MASSLOSS == "FIX" ? @sprintf("%4.2f", f_def) : ""

    #^^^ MODEL PARAMETERS
    ####
    #VVV RUN



    f_B_system = get_multiplicity(fB; mi=0, mf=100., custom_value=fB_value)
    f_B(x) = f_B_system(x)   #f_B_system(x) ./ (2. - f_B_system(x))
    f_S(x) = 1-f_B(x)

    integral_logP, _ = quadgk(logp -> SANA_logP_pdf(logp, pow_p), minimum(logPs)-delta, maximum(logPs)+delta, rtol=1e-8)
    integral_q, _    = quadgk(q    -> SANA_q_pdf(q, pow_q),    minimum(qs)-delta, maximum(qs)+delta, rtol=1e-9)
    integral_logM, _ = integrate_logM(minimum(logMs)-delta, maximum(logMs)+delta,  pow_m)

    function block_logP_pdf(logP)
        integral_block, _ = quadgk(x -> SANA_logP_pdf(x, pow_p), logP-delta, logP+delta, rtol=1e-8)
        return integral_block/integral_logP
    end
    function block_q_pdf(q)
        integral_block, _ = quadgk(x -> SANA_q_pdf(x, pow_q), q-delta, q+delta, rtol=1e-9)
        return integral_block/integral_q
    end
    function block_logM_pdf(logM)
        integral_block, _ = integrate_logM(logM-delta, logM+delta, pow_m)
        return integral_block/integral_logM
    end
    function block_logM_BINARY_pdf(logM)
        integral_block, _ = integrate_logM_convolve(logM-delta, logM+delta, pow_m, f_B)
        return integral_block/integral_logM
    end
    function block_logM_SINGLE_pdf(logM)
        integral_block, _ = integrate_logM_convolve(logM-delta, logM+delta, pow_m, f_S)
        return integral_block/integral_logM
    end

    
    pdfs = SortedDict("logM"=>Dict(), "logP" => Dict(), "q"=>Dict())
    for logM in logMs
        logM_t = to_key(logM) 
        pdfs["logM"][logM_t] =  block_logM_BINARY_pdf(logM)
    end
    for q in qs
        q_t = to_key(q) 
        pdfs["q"][q_t] = block_q_pdf(q)
    end
    for logP in logPs
        logP_t = to_key(logP) 
        pdfs["logP"][logP_t] =  block_logP_pdf(logP)
    end

    SINGLES = Dict()
    for logM in logMs
        logM_t = to_key(logM) 
        SINGLES[logM_t]=Dict()
        SINGLES[logM_t]["pdf"]=0.
        SINGLES[logM_t]["SN"]=Dict(:WD =>0., :IIP =>0., :IIb =>0., :SN87A=>0., 
                                    :Ibc=>0., :BH =>0., :IIn =>0., :Ibn =>0., :X=>0.)
    end
    for modkey in keys(MODELS)
        logM_t = get_from_modkey(what="logM", key=modkey, returntype=String)
        logP_t = get_from_modkey(what="logP", key=modkey, returntype=String)
        q_t    = get_from_modkey(what="q",    key=modkey, returntype=String)
        pdfs_p = pdfs["logP"][logP_t]
        pdfs_q = pdfs["q"][q_t]
        pdfs_m = pdfs["logM"][logM_t]
        MODELS[modkey]["pdf"] = (q     = pdfs_q,
                                 logP  = pdfs_p,
                                 logM  = pdfs_m, 
                                 logP_q= pdfs_p*pdfs_q,
                                 tot   = pdfs_p*pdfs_q*pdfs_m)
    end





    #sdelta=0.020/2
    sdelta = 0.0001
    for logm in minimum(logMs)-delta:sdelta:maximum(logMs)+delta
        logm_t = @sprintf("%.5f", logm)
        m = 10 ^ logm
        logm_min = max(minimum(logMs)-delta, logm-sdelta/2)
        logm_max = min(maximum(logMs)+delta, logm+sdelta/2)

        sn=Dict(:WD =>0., :IIP =>0., :IIb =>0., :SN87A=>0., 
        :Ibc=>0., :BH =>0., :IIn =>0., :Ibn =>0., :X=>0.)

        sn_type, outcome = SN_OUTPUT_manual(interpol_single_Mmax["M_end"](m),
                            interpol_single_Mmax["M_hecore_end"](m), 
                            interpol_single_Mmax["M_cocore_end"](m), 
                            m, 
                            interpol_single_Mmax["M_env_end"](m),
                            interpol_single_Mmax["X_C"](m),
                            "end",
                            EXP_CRIT)
        sn[sn_type] = 1.
        counter = 0
        for logM in logMs
            ! (logm_min < logM+delta) && continue 
            ! (logm_max > logM-delta) && continue 
            logm_low  = max(logm_min, logM-delta)
            logm_high = min(logm_max, logM+delta)

            logM_t = @sprintf("%.2f", logM)
            counter += 1
            (counter > 2) && (throw(ErrorException("SUB-BLOCKS CONTRIBUTING TO MORE THAN ONE BLOCK!")))  

            pdf = integrate_logM_convolve(logm_low, logm_high, pow_m, f_S)[1]/integral_logM
            SINGLES[logM_t]["pdf"] += pdf
            for outcome in outcomes 
                SINGLES[logM_t]["SN"][outcome] += sn[outcome] * pdf
            end 
        end 

    end
    for logM in logMs
        logM_t = @sprintf("%.2f", logM)
        for outcome in outcomes 
            if fB_value >= 1
                SINGLES[logM_t]["SN"][outcome] = 0
                SINGLES[logM_t]["pdf"] = 0
            else
                SINGLES[logM_t]["SN"][outcome] /= SINGLES[logM_t]["pdf"]

            end

        end 
    end 



    # integral_manual = 0
    # for logP in range(minimum(logPs), maximum(logPs), step=0.05)
    #     block = block_logP_pdf(logP)
    #     integral_manual += block
    #     @printf("logP=%4.2f - p = %5.1f%%,  cumul = %6.1f%%\n", logP, block*100, 100*integral_manual)
    # end
    # @printf("logP blocks total pdf: %.3f\n", integral_manual)

    stats = Dict()
    stats["1"] = Dict()
    stats["2"] = Dict()
    stats["M"] = Dict()
    stats["S"] = Dict()

    for outcome in outcomes
        stats["1"][outcome] = 0
        stats["2"][outcome] = 0
        stats["M"][outcome] = 0
        stats["S"][outcome] = 0
    end

    abs_sn_number_zams = 0
    abs_number_zams = 0
    abs_sn_number_zams_s = 0
    abs_number_zams_s = 0

    for modkey in keys(MODELS)
        mod = MODELS[modkey]
        pdf =  mod["pdf"].tot
        
        firstSN = mod["1_whichstar"]
        primary   = (firstSN == "1") ? "1" : "2"
        secondary = (firstSN == "1") ? "2" : "1"

        abs_number_zams += 2*pdf
        for outcome in outcomes
            post_sn_runs1 = length(mod["1"].SN)
            post_sn_runs2 = length(mod["2"].SN)
            sn1 = sum( [ (mod["1"].SN[s] == outcome ? 1 : 0) * 1                                                                       for s in range(1,post_sn_runs1)])/post_sn_runs1
            sn2 = sum( [ (mod["2"].SN[s] == outcome ? 1 : 0) * (s == mod["unboundstate"]["ix"] ? mod["unboundstate"]["counter"] : 1 )  for s in range(1,post_sn_runs2)])/( post_sn_runs2 + max(mod["unboundstate"]["counter"],1)-1 )
            if mod["1_whichstar"] == "M"
                stats["M"][outcome] += sn1 * pdf
                stats["1"][:X]     += sn1 * pdf
                stats["2"][:X]     += sn2 * pdf
                (outcome in SN_outcomes && sn1>0.) && (abs_sn_number_zams  += sn1 * pdf)
            else
                stats["M"][:X]     += sn1 * pdf
                stats["1"][outcome] += (primary == firstSN ? sn1 : sn2) * pdf
                stats["2"][outcome] += (primary == firstSN ? sn2 : sn1) * pdf
                (outcome in SN_outcomes && sn1>0.) && (abs_sn_number_zams  += sn1 * pdf)
                (outcome in SN_outcomes && sn2>0.) && (abs_sn_number_zams  += sn2 * pdf)
            end 
        end
    end
    for logM in logMs
        logM_t = to_key(logM)
        abs_number_zams_s += SINGLES[logM_t]["pdf"]
        for outcome in outcomes
            stats["S"][outcome] += SINGLES[logM_t]["SN"][outcome] * SINGLES[logM_t]["pdf"]
            (outcome in SN_outcomes && SINGLES[logM_t]["SN"][outcome]>0. ) && (abs_sn_number_zams_s  += SINGLES[logM_t]["SN"][outcome] * SINGLES[logM_t]["pdf"])
        end
    end


    check = Dict()
    for prog in ["1", "2", "M"]
        check[prog] = 0
        check["S"] = 0
        for outcome in outcomes
            check[prog] += stats[prog][outcome]
            check["S"] += stats["S"][outcome]
        end
        
        if ! (abs(check[prog]+check["S"]-1)< 0.001)
            for outcome in outcomes
                @printf("[%2s] - %5s  :  %7.4f + %7.4f = %7.4f\n", prog, outcome, stats[prog][outcome], stats["S"][outcome],  stats[prog][outcome]+stats["S"][outcome])
            end
            @printf("%7.4f (%2s) + %7.4f (%2s) = %7.4f\n", check[prog], prog, check["S"], "S", check[prog]+check["S"])
            throw(ErrorException("stats are not unitary"))
        end
    end

    total_SNe = 0
    for out in SN_outcomes 
        total_SNe += stats["1"][out]+stats["2"][out]+stats["M"][out]+stats["S"][out]
    end

    stats_SN = Dict()
    stats_BH = Dict()
    for prog in progenitors
        stats_SN[prog] = Dict()
        stats_BH[prog] = Dict()

    end 
    for SN_outcome in SN_outcomes 
        stats_SN["1"][SN_outcome] = stats["1"][SN_outcome]/total_SNe
        stats_SN["2"][SN_outcome] = stats["2"][SN_outcome]/total_SNe
        stats_SN["M"][SN_outcome] = stats["M"][SN_outcome]/total_SNe
        stats_SN["S"][SN_outcome] = stats["S"][SN_outcome]/total_SNe
    end
    stats_BH["1"] = stats["1"][:BH]/total_SNe
    stats_BH["2"] = stats["2"][:BH]/total_SNe
    stats_BH["M"] = stats["M"][:BH]/total_SNe
    stats_BH["S"] = stats["S"][:BH]/total_SNe

    total_SNe2 = 0
    stats_SN_tot = Dict()
    for SN_outcome in SN_outcomes 
        stats_SN_tot[SN_outcome] = sum( [stats_SN[prg][SN_outcome] for prg in progenitors] )
        total_SNe2 += stats_SN_tot[SN_outcome]
    end 
    stats_BH_tot = sum( [stats_BH[prg] for prg in progenitors] )

    ! (0.9999 < total_SNe2 < 1.0001) ? throw(ErrorException("SNe not normalized")) : nothing








    @printf("%10s | %7s %7s %7s %7s %7s %7s (%7s)\n", "which",  "IIP", "SN87A", "IIb", "Ibc", "IIn", "Ibn" , "BH")
    for prg in progenitors
        stats_SN[prg]["total"] = sum([stats_SN[prg][SN_outcome] for SN_outcome in SN_outcomes])
        @printf("%10s | %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f (%7.2f) | %7.2f\n", prg,  
                100*stats_SN[prg][:IIP], 100*stats_SN[prg][:SN87A], 100*stats_SN[prg][:IIb],
                100*stats_SN[prg][:Ibc], 100*stats_SN[prg][:IIn],   100*stats_SN[prg][:Ibn],  100*stats_BH[prg], 100*stats_SN[prg]["total"])

    end 
    @printf("%10s | %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f (%7.2f)\n", "TOTAL",  
               100*stats_SN_tot[:IIP], 100*stats_SN_tot[:SN87A], 100*stats_SN_tot[:IIb],
               100*stats_SN_tot[:Ibc], 100*stats_SN_tot[:IIn],  100*stats_SN_tot[:Ibn], 100*stats_BH_tot)



    if WRITE_TO_FILE
        tot_SN_CORR = out_dir*"tot_SN_CORR.dat" 
        isfile(tot_SN_CORR) ? i=1 : i=0 
        fp = open(tot_SN_CORR, "a")
        if i ==0 
            write(fp,  @sprintf("%10s ", "MER_CRIT"))
            write(fp,  @sprintf("%10s ", "MER_dM"))
            write(fp,  @sprintf("%10s ", "MER_EOL"))
            write(fp,  @sprintf("%10s ", "EXP_CRIT"))
            write(fp,  @sprintf("%10s ", "fB"))
            write(fp,  @sprintf("%10s ", "pow_logM"))
            write(fp,  @sprintf("%10s ", "pow_logP"))
            write(fp,  @sprintf("%10s ", "pow_q"))
            write(fp,  @sprintf("%10s ", "total_SNe"))
            write(fp,  @sprintf("%10s ", "zams_tot"))
            write(fp,  @sprintf("%10s ", "sn_tot"))
            write(fp,  @sprintf("%10s ", "zams_b"))
            write(fp,  @sprintf("%10s ", "sn_b"))
            write(fp,  @sprintf("%10s ", "zams_s"))
            write(fp,  @sprintf("%10s ", "sn_s"))
            write(fp,  @sprintf("%10s ", "IIP"))
            [write(fp,  @sprintf("%10s ", "IIP_"*prg)) for prg in progenitors ] 
            write(fp,  @sprintf("%10s ", "SN87A"))
            [write(fp,  @sprintf("%10s ", "SN87A_"*prg)) for prg in progenitors ] 
            write(fp,  @sprintf("%10s ", "IIb"))
            [write(fp,  @sprintf("%10s ", "IIb_"*prg)) for prg in progenitors ] 
            write(fp,  @sprintf("%10s ", "Ibc"))
            [write(fp,  @sprintf("%10s ", "Ibc_"*prg)) for prg in progenitors ] 
            write(fp,  @sprintf("%10s ", "IIn"))
            [write(fp,  @sprintf("%10s ", "IIn_"*prg)) for prg in progenitors ] 
            write(fp,  @sprintf("%10s ", "Ibn"))
            [write(fp,  @sprintf("%10s ", "Ibn_"*prg)) for prg in progenitors ] 
            write(fp,  @sprintf("%10s ", "BH"))
            [write(fp,  @sprintf("%10s ", "BH_"*prg)) for prg in progenitors ] 
            write(fp,  @sprintf("\n"))
            i=1
        end 
        write(fp,  @sprintf("%10s ", MERGER_CRITERION))
        write(fp,  @sprintf("%10s ", MERGER_MASSLOSS))
        write(fp,  @sprintf("%10s ", MERGER_EOL))
        write(fp,  @sprintf("%10s ", EXP_CRIT))
        write(fp,  fB == :interpolate ? @sprintf("%10s ", fB_value) : @sprintf("%10.3f ", fB_value))
        write(fp,  @sprintf("%10.3f ", pow_m))
        write(fp,  @sprintf("%10.3f ", pow_p))
        write(fp,  @sprintf("%10.3f ", pow_q))
        write(fp,  @sprintf("%10.4f ", total_SNe))
        write(fp,  @sprintf("%10.4f ", abs_number_zams+abs_number_zams_s))
        write(fp,  @sprintf("%10.4f ", abs_sn_number_zams+abs_sn_number_zams_s))
        write(fp,  @sprintf("%10.4f ", abs_number_zams))
        write(fp,  @sprintf("%10.4f ", abs_sn_number_zams))
        write(fp,  @sprintf("%10.4f ", abs_number_zams_s))
        write(fp,  @sprintf("%10.4f ", abs_sn_number_zams_s))
        write(fp,  @sprintf("%10.3f ",  100 *stats_SN_tot[:IIP]))
        [write(fp, @sprintf("%10.3f ",  100*stats_SN[prg][:IIP])) for prg in progenitors ] 
        write(fp,  @sprintf("%10.3f ",  100 *stats_SN_tot[:SN87A]))
        [write(fp, @sprintf("%10.3f ",  100*stats_SN[prg][:SN87A])) for prg in progenitors ] 
        write(fp,  @sprintf("%10.3f ",  100 *stats_SN_tot[:IIb]))
        [write(fp, @sprintf("%10.3f ",  100*stats_SN[prg][:IIb])) for prg in progenitors ] 
        write(fp,  @sprintf("%10.3f ",  100 *stats_SN_tot[:Ibc]))
        [write(fp, @sprintf("%10.3f ",  100*stats_SN[prg][:Ibc])) for prg in progenitors ] 
        write(fp,  @sprintf("%10.3f ",  100 *stats_SN_tot[:IIn]))
        [write(fp, @sprintf("%10.3f ",  100*stats_SN[prg][:IIn])) for prg in progenitors ] 
        write(fp,  @sprintf("%10.3f ",  100 *stats_SN_tot[:Ibn]))
        [write(fp, @sprintf("%10.3f ",  100*stats_SN[prg][:Ibn])) for prg in progenitors ] 
        write(fp,  @sprintf("%10.3f ",  100 *stats_BH_tot))
        [write(fp, @sprintf("%10.3f ",  100*stats_BH[prg])) for prg in progenitors ] 
        write(fp,  @sprintf("\n"))

        close(fp)
    end

    if only_SN 
        @printf("Removing models that do not produce any SN... n. models = %6d -> ", length(keys(MODELS)))
        for mod in keys(MODELS)
            if MODELS[mod]["1"].SN[1] ∈ [:WD, :X]  && MODELS[mod]["2"].SN[1] ∈ [:X, :WD] && (length(MODELS[mod]["1"].SN) + length(MODELS[mod]["2"].SN) == 2)
                # @printf("Removing %s from MODELS since it does not produce any SN\n", mod)
                delete!(MODELS, mod)
            end
        end
        @printf(" %6d \n", length(keys(MODELS)))

    end

    MODELS, SN_DATA, BH_binaries = reorganize_SN_data(;
                                            file=f,
                                            MODELS=MODELS, 
                                            f_B=fB_value,
                                            EXP_CRIT = EXP_CRIT,
                                            pow_m = pow_m,
                                            total_SNe  =  total_SNe,
                                            only_SN = only_SN,
                                            logMs = logMs)


    return MODELS, SN_DATA, BH_binaries, total_SNe
end
                



function reorganize_SN_data(;   file = nothing,
                                MODELS=nothing, 
                                f_B=   nothing, 
                                EXP_CRIT  = nothing, 
                                pow_m     = nothing, 
                                total_SNe = nothing,
                                only_SN   = false,
                                logMs = 0.07:0.05:2.00)


    for model in keys(MODELS)
        mod=MODELS[model]
        for s in range(1,length(mod["1"].SN))
            if mod["1"].SN[s] == :Ibn 
                if mod["1"].dM_C[s] > 10 
                    mod["1"].dM_C[s] = (file.M_he_core_1_end_BURN_He_1[mod["mod"]] - mod["1"].endvals[s].M_end) - 0.1
                end
            end
        end
    end
    


    f_single = CSV.read(SG_dir*"CORR_GRID_single_150.dat", 
                        header=1, DataFrame,delim=' ',ignorerepeated=true)
    deleteat!(f_single, findall( !>(0) , f_single.M_he_core_end_BURN_He) )


    f_S = 1-f_B
    delta =0.050/2
    integral_logM, _ = integrate_logM( minimum(logMs)-delta, maximum(logMs)+delta, pow_m )
    dump = 0

    sdelta = 0.0001
    SINGLES = SortedDict()
    for logm in minimum(logMs)-delta:sdelta:maximum(logMs)+delta
        f_S <= 0 && continue
        logm_t = @sprintf("%.5f", logm)
        SINGLES[logm_t] = Dict()

        m = 10 ^ logm
        logm_min = max(minimum(logMs)-delta, logm-sdelta/2)
        logm_max = min(maximum(logMs)+delta, logm+sdelta/2)
        pdf = integrate_logM(logm_min, logm_max, pow_m)[1] / integral_logM * f_S 

        SINGLES[logm_t] = Dict()
        SINGLES[logm_t]["pdf"]  = (tot=pdf,)
        SINGLES[logm_t]["mod"] = NaN
        SINGLES[logm_t]["casus"] = "SINGLE"
        SINGLES[logm_t]["ECI"] = [ECI.estimate_ECI(Mej=NaN, Eexp=NaN, a=+Inf, M2=NaN, R2=NaN, RL2=+Inf )]
        SINGLES[logm_t]["post_SN_orbit"] = [(a=Inf, a_peri=Inf, P=Inf, e=Inf, ϵ=Inf, 
                                            orbit = :unbound, 
                                            will_one_periastron_occur=false, 
                                            when_will_periastron_occur = NaN)]
        sn, endvals = SN_OUTPUT_manual(interpol_single_Mmax["M_end"](m),
                                    interpol_single_Mmax["M_hecore_end"](m), 
                                    interpol_single_Mmax["M_cocore_end"](m), 
                                    m,
                                    interpol_single_Mmax["M_env_end"](m),
                                    interpol_single_Mmax["X_C"](m),
                                    "end",
                                    EXP_CRIT)
        counter = 0
        SINGLES[logm_t]["1"] = (SN      = [sn],
                                Mmax    = [m],
                                Mend    = [interpol_single_Mmax["M_end"](m)],
                                Mhe     = [interpol_single_Mmax["M_hecore_end"](m)],
                                Mco     = [interpol_single_Mmax["M_cocore_end"](m)],
                                Xc      = [interpol_single_Mmax["X_C"](m)],
                                dM_C    = [NaN],
                                endvals = [endvals],
                                t_SN    = [f_age_CB(m)],
                                case    = [""]
                                )
        SINGLES[logm_t]["1_whichstar"]  = "S"
        SINGLES[logm_t]["pre_SN_orbit"] = (a=NaN, v1=NaN, v2=NaN)
        SINGLES[logm_t]["1st_kick"] = Float32[NaN]
        SINGLES[logm_t]["1st_kick_whole"] = [(NaN,NaN,NaN)]


    end

    counter = 0 
    counterM = 0 
    counterP = 0 
    counterq = 0 

    quantities_to_collect  = ["t", "Mf", "MM", "pdf", "Mi", "E_exp", 
                            "M_rem", "M_ni", "Type", "dM_C", "M_rem_b", "orbSN1",
                            "model_n", "case", "Mej", "Mhe", "Menv", "Mheshell", "Mco", "ECI_ogata", "ECI_manual", "ECI_tau",
                            "postkick_orbit", "periastron_occurs?",  "when_periastron", "postkick_period", "firstSN"]
    progenitors_to_collect = ["tot", "1", "2", "M", "S", "1-MT", "2-MT", "1-NI", "2-NI", "S+NI"]
    SN_DATA =     Dict()
    for qt in quantities_to_collect
        SN_DATA[qt] = Dict()
        for prg in progenitors_to_collect 
            SN_DATA[qt][prg] = qt in ["Type", "case", "ECI_ogata", "ECI_manual", "postkick_orbit", "periastron_occurs?", "firstSN"] ? [] : Float32[] 
        end
    end
    SN_DATA["tot_SNe"] = total_SNe

    binary_BHs_n = []

    counter = 0
    counter2 = 0

    function pusher(SN_DATA,  prg, whichSN, mod, ix,  logM, q)

        # sum( [ (mod["2"].SN[s] == outcome ? 1 : 0) * (s == mod["unboundstate"]["ix"] ? mod["unboundstate"]["counter"] : 1 )  for s in range(1,post_sn_runs2)])/( post_sn_runs2+ max(mod["unboundstate"]["counter"]-1,1) )
        unboundstates = (whichSN=="2") ? max(mod["unboundstate"]["counter"],1) : 1
        n_sample = length(mod[whichSN].SN)+unboundstates-1
        unboundstateboost = (whichSN=="2" && ix == mod["unboundstate"]["ix"]) ? mod["unboundstate"]["counter"] : 1
        mod_n = mod["mod"]
        this_mod = mod[whichSN]
        push!( SN_DATA["t"][prg],        this_mod.t_SN[ix])
        push!( SN_DATA["model_n"][prg],  mod_n)
        push!( SN_DATA["case"][prg],     this_mod.case[ix])
        push!( SN_DATA["Type"][prg],     this_mod.SN[ix])
        push!( SN_DATA["Mf"][prg],       this_mod.endvals[ix].M_end )
        push!( SN_DATA["Mhe"][prg],      this_mod.endvals[ix].M_he )
        push!( SN_DATA["Mco"][prg],      this_mod.endvals[ix].M_co )
        push!( SN_DATA["Menv"][prg],     this_mod.endvals[ix].M_end - this_mod.endvals[ix].M_he )
        push!( SN_DATA["Mheshell"][prg], this_mod.endvals[ix].M_he  - this_mod.endvals[ix].M_co )
        push!( SN_DATA["MM"][prg],       this_mod.Mmax[ix] )
        push!( SN_DATA["dM_C"][prg],     this_mod.dM_C[ix])
        push!( SN_DATA["pdf"][prg],      mod["pdf"].tot/total_SNe/n_sample * unboundstateboost)
        push!( SN_DATA["E_exp"][prg],    this_mod.endvals[ix].E_exp/1e51)
        push!( SN_DATA["M_rem"][prg],    this_mod.endvals[ix].M_remnant_g )
        push!( SN_DATA["M_rem_b"][prg],  this_mod.endvals[ix].M_remnant_b )
        push!( SN_DATA["M_ni"][prg],     this_mod.endvals[ix].M_ni )
        push!( SN_DATA["Mej"][prg],      this_mod.endvals[ix].M_end -this_mod.endvals[ix].M_remnant_b )
        push!( SN_DATA["firstSN"][prg],  whichSN ==  mod["1_whichstar"])
        push!( SN_DATA["ECI_ogata"][prg],          mod["ECI"][ix]["ogata"]["ECI?"] )
        push!( SN_DATA["ECI_tau"][prg],            mod["ECI"][ix]["ogata"]["tau_ECI"] )
        push!( SN_DATA["ECI_manual"][prg],         mod["ECI"][ix]["manual"]["ECI?"])
        push!( SN_DATA["postkick_orbit"][prg],     mod["post_SN_orbit"][ix].orbit)
        push!( SN_DATA["periastron_occurs?"][prg], mod["post_SN_orbit"][ix].will_one_periastron_occur)
        push!( SN_DATA["when_periastron"][prg],    mod["post_SN_orbit"][ix].when_will_periastron_occur/yr)
        push!( SN_DATA["postkick_period"][prg],    mod["post_SN_orbit"][ix].P/day)

        push!( SN_DATA["Mi"][prg], 10. ^ logM * (prg == "2" ? q : 1))
    end 


    for modkey in keys(MODELS)
        mod = MODELS[modkey]
        logM = get_from_modkey(what="logM", key=modkey)
        q = get_from_modkey(what="q", key=modkey)

        if mod["1_whichstar"] == "M"
            if mod["1"].SN[1] in [:WD, :X]
                one = sum( [ (mod["2"].SN[s] == :X ? 1 : 0)  for s in range(1,length(mod["2"].SN))] )
                one != 1 && @printf("sum SN2[X] = %.3f : diff1 = %.3e\n", one, abs(1-one))
                one != 1 && throw(ErrorException("Check merger: got 1=M but 2 is not X!"))
                continue
            end
            pusher(SN_DATA, "M",   "1", mod, 1, logM, q)
            pusher(SN_DATA, "tot", "1", mod, 1, logM, q)
        else
            firstSN = mod["1_whichstar"]
            primary   = (firstSN == "1") ? "1" : "2"
            secondary = (firstSN == "1") ? "2" : "1"

            for s in range(1,length(mod[primary].SN))
                RLOF_pre1stSN = length(mod["1"].case[1]) > 0 
                RLOF =  if (primary == firstSN) 
                            RLOF_pre1stSN 
                        else
                            RLOF_pre1stSN || length(mod["2"].case[s])>0
                        end
                # if ! (mod[primary].SN[s] in [:X, :WD] )
                    pusher(SN_DATA, "1",   primary, mod, s, logM, q)
                    pusher(SN_DATA, "tot", primary, mod, s, logM, q)
                    if RLOF
                        pusher(SN_DATA, "1-MT", primary, mod, s, logM, q)
                    else
                        pusher(SN_DATA, "1-NI", primary, mod, s, logM, q)
                        pusher(SN_DATA, "S+NI", primary, mod, s, logM, q)
                    end
                # end 
            end

            for s in range(1,length(mod[secondary].SN))
                RLOF_pre1stSN = length(mod["1"].case[1]) > 0 
                RLOF =  if (secondary == firstSN) 
                            RLOF_pre1stSN 
                        else
                            RLOF_pre1stSN || length(mod["2"].case[s])>0
                        end
                # if ! (mod[secondary].SN[s] in [:X, :WD] )
                    pusher(SN_DATA, "2",   secondary, mod, s, logM, q)
                    pusher(SN_DATA, "tot", secondary, mod, s, logM, q)
                    if RLOF
                        pusher(SN_DATA, "2-MT", secondary, mod, s, logM, q)
                    else
                        pusher(SN_DATA, "2-NI", secondary, mod, s, logM, q)
                        pusher(SN_DATA, "S+NI", secondary, mod, s, logM, q)
                    end
                # end 
            end

            for p in range(1,length(mod[primary].SN))
                for s in range(1,length(mod[secondary].SN))
                    ((mod[primary].SN[p] == :BH )&& (mod[secondary].SN[s] == :BH)) && push!(binary_BHs_n, mod)
                end
            end
        end
    end


    for logm in keys(SINGLES)
        (only_SN && SINGLES[logm]["1"].SN[1] ∈ [:X, :WD]) && continue
        MODELS[logm*"_s"]=SINGLES[logm]
        mod = SINGLES[logm]
        logM = parse(Float64 ,logm)

        if isnan(mod["1"].endvals[1].M_end)
            counter2 += 1 
            throw(ErrorException("why is this nan?"))
        end

        pusher(SN_DATA, "S",    "1", mod, 1, logM, 1)
        pusher(SN_DATA, "S+NI", "1", mod, 1, logM, 1)
        pusher(SN_DATA, "tot",  "1", mod, 1, logM, 1)

    end 

    function gix(prg::String, type)
        type isa Symbol && (type = [type])
        return [ix for (ix, t) in enumerate(SN_DATA["Type"][prg]) if (type == [:all] || t in type)]
    end 


    # total_SNe2 = sum(SN_DATA["pdf"]["tot"])*100
    # ! (0.9999 < total_SNe2 < 1.0001) ? throw(@printf(ErrorException("SNe not normalized")) : nothing
    @printf("prg %6s %6s %6s %6s %6s %6s :%6s |%6s\n", "IIP", "SN87A", "IIb", "Ibc", "IIn", "Ibn", "BH", "Total")
    for prg in ["1", "2", "M", "S"]
        @printf("%4s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f :%6.2f | %6.2f\n", prg,
        sum(SN_DATA["pdf"][prg][gix(prg,:IIP)])*100,  sum(SN_DATA["pdf"][prg][gix(prg,:SN87A)])*100, 
        sum(SN_DATA["pdf"][prg][gix(prg,:IIb)])*100,  sum(SN_DATA["pdf"][prg][gix(prg,:Ibc)])*100,
        sum(SN_DATA["pdf"][prg][gix(prg,:IIn)])*100, sum(SN_DATA["pdf"][prg][gix(prg,:Ibn)])*100, 
        sum(SN_DATA["pdf"][prg][gix(prg,:BH)])*100, 100*sum(SN_DATA["pdf"][prg][gix(prg,[:IIP, :SN87A, :IIb, :Ibc, :IIn, :Ibn])]) )
    end
    @printf("Placing non-interacting systems inside single-stars.\n")
    for prg in ["1-MT", "2-MT", "M", "S+NI"]
        @printf("%4s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f :%6.2f | %6.2f\n", prg,
        sum(SN_DATA["pdf"][prg][gix(prg,:IIP)])*100,  sum(SN_DATA["pdf"][prg][gix(prg,:SN87A)])*100, 
        sum(SN_DATA["pdf"][prg][gix(prg,:IIb)])*100,  sum(SN_DATA["pdf"][prg][gix(prg,:Ibc)])*100,
        sum(SN_DATA["pdf"][prg][gix(prg,:IIn)])*100, sum(SN_DATA["pdf"][prg][gix(prg,:Ibn)])*100, 
        sum(SN_DATA["pdf"][prg][gix(prg,:BH)])*100, 100*sum(SN_DATA["pdf"][prg][gix(prg,[:IIP, :SN87A, :IIb, :Ibc, :IIn, :Ibn])]) )
    end


    
    GC.gc() 

    return MODELS, SN_DATA, binary_BHs_n



end




