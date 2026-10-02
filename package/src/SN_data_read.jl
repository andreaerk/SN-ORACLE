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



function get_idx(prg::String, type; database=SN_DATA, filter = nothing)
    out = Int[]
    type isa Symbol && (type = [type])
    (filter === nothing) && (filter = [true for i in range(1,length(database["model_n"][prg]))])
    (length(filter) != length(database["model_n"][prg])) && throw(ErrorException("[get_idx] filter not of the same length as database array!"))
    for ix in range(1,length(database["model_n"][prg]))
        !(type == [:all] || database["Type"][prg][ix] in type) && continue
        !filter[ix] && continue
        push!( out, ix )
    end
    
    return out
end 

function remove_values(arr, vals)
    vals_set = Set(vals)
    return [x for x in arr if !(x in vals_set)]
end
function remove_at_indices(indices, arrays...)
    return ( [arr[i] for i in eachindex(arr) if !(i in indices)] for arr in arrays )
end
function get_data_sne(database, prg, type, what; filter=nothing)
    (filter === nothing) && (filter = [true for i in range(1,length(database["model_n"][prg]))])
    return database[what][prg][get_idx(prg, type; database=database, filter=filter)]
end
function remove_values_and_linked(vals_to_remove, arr, arrays...)
    indices_to_remove = findall(x -> x in Set(vals_to_remove), arr)
    return ( [a[i] for i in eachindex(arr) if !(i in indices_to_remove)] for a in (arr, arrays...) )
end
function remove_values_and_linked_new(vals_to_remove, arr, arrays...)
    S = Set(vals_to_remove)
    mask = .!(in.(arr, Ref(S)))   # Boolean mask: true = keep
    return (a[mask] for a in (arr, arrays...))
end


function report_stats(vals, pdf; label = "", print2screen = false, newline = false, unit = (1, ""), format_str = "%4.1f", n_peaks = 1, multip1_threshold=[15], multip_percentiles=0:5:100)
    pdfs = pdf
    report_stats_common(vals, pdfs, label; print2screen, newline, unit, format_str, n_peaks, multip1_threshold, multip_percentiles)
end

function report_stats(database, what, prg, ttp; filter = nothing, print2screen = false, newline = false, unit = (1, ""), format_str = "%4.1f", n_peaks = 1, multip1_threshold=[15], multip_percentiles=0:5:100)
    val = nothing 
    pdfs = nothing 
    if filter === nothing
        vals =  [database[what][prg][ix]/unit[1] for ix in get_idx(prg, ttp; database = database)]
        pdfs = [database["pdf"][prg][ix] for ix in get_idx(prg, ttp; database = database)]
    else
        vals =  [database[what][prg][ix]/unit[1] for ix in get_idx(prg, ttp; database = database, filter=filter)]
        pdfs = [database["pdf"][prg][ix] for ix in get_idx(prg, ttp; database = database, filter=filter)]
    end

    report_stats_common(vals, pdfs, what; print2screen, newline, unit, format_str, n_peaks, multip1_threshold, multip_percentiles)
end

function report_stats_common(vals, pdfs, label; print2screen, newline, unit, format_str, n_peaks, multip1_threshold, multip_percentiles)
    vals, pdfs = remove_values_and_linked_new([NaN, Inf, -Inf], vals, pdfs)
    pdfs = Weights(pdfs)
    val_min = minimum(vals) 
    val_max = maximum(vals)
    output_vals = []

    avg, std = mean_and_std( vals, pdfs )

    q_05, q_95, q_50 = quantile(vals, pdfs, 0.05), quantile(vals, pdfs, 0.95),  quantile(vals, pdfs, 0.50)

    print2screen && (fmt_string = @sprintf("%s : %s<%s[+%s-%s]<%s [%s]. %s", "%5s", format_str, format_str, format_str, format_str, format_str, "%5s", "%s"))
    print2screen && (fmt = Printf.Format(fmt_string))

    print2screen && Printf.format(stdout, fmt, label, val_min/unit[1], avg/unit[1],  (q_95-avg)/unit[1], (avg-q_05)/unit[1], val_max/unit[1], unit[2], (newline ? "\n" : "") )
    
    print2screen && @printf("\nPerc. Analysis\n")
    n=0
    if print2screen 

        for qt in multip_percentiles
            n_digs = format_str[2] 
            fmt_string = @sprintf("%s%s%s ", "%", n_digs, "s")
            fmt = Printf.Format(fmt_string)
            print2screen &&  Printf.format(stdout, fmt, qt)
            n+= 1
            (n%15 == 0) && @printf("\n")
        end
        print2screen && @printf("\n")

        n=0
        for qt in multip_percentiles
            fmt = Printf.Format(format_str*" ")
            print2screen &&  Printf.format(stdout, fmt, quantile(vals,pdfs,qt/100)/unit[1])
            n+= 1
            (n%15 == 0) && @printf("\n")
        end
        print2screen && @printf("\n")
        @printf("\n")
    end 

    
    if n_peaks > 1 
         
        print2screen && @printf("\n")
        qtl_thresholds = copy(multip1_threshold)
        push!(qtl_thresholds, 100)
        print2screen  && @printf("Having chosen %s quantile as the threshold for the multipeak, here are the statistics\n", multip1_threshold)
        # peaks = [[vals[vals .<  thershold], pdfs[vals  .<  thershold] ]  ,
        #          [vals[vals .>= thershold], pdfs[vals  .>= thershold] ]  ]
        threshold_old = 0
        for quantile_threshold in qtl_thresholds 
            
            thershold = quantile(vals,pdfs, quantile_threshold/100)

            vals_pk = vals[threshold_old .< vals .<=  thershold]
            pdfs_pk = pdfs[threshold_old .< vals .<=  thershold]

            avg, std = mean_and_std( vals_pk, pdfs_pk )
            val_min = minimum(vals_pk) 
            val_max = maximum(vals_pk)
            q_05, q_95 = quantile(vals_pk, pdfs_pk, 0.05), quantile(vals_pk, pdfs_pk, 0.95)
           
            print2screen && (fmt_string = @sprintf("%s - %s<%s[+%s-%s]<%s [%s]- %s", "%5s", format_str, format_str, format_str, format_str, format_str, "%5s", "%s"))
            print2screen && (fmt = Printf.Format(fmt_string))
            print2screen && Printf.format(stdout, fmt, what, val_min/unit[1], avg/unit[1],  (q_95-avg)//unit[1], (avg-q_05)//unit[1], val_max//unit[1], unit[2],  "\n"  )
        
            push!(output_vals, [avg, q_95-avg, avg-q_05, q_50])

            threshold_old = thershold
        end


    else 

        push!(output_vals, [avg, q_95-avg, avg-q_05, q_50])

    end

    return output_vals 


end

function weighted_quantile_gap_method(database, what, prg, ttp; filter = nothing, unit = (1, ""))
    val = nothing 
    pdfs = nothing 
    if filter === nothing
        vals =  [database[what][prg][ix]/unit[1] for ix in get_idx(prg, ttp)]
        pdfs = Weights([database["pdf"][prg][ix] for ix in get_idx(prg, ttp)])
    else
        vals =  [database[what][prg][ix]/unit[1] for ix in get_idx(prg, ttp, filter=filter)]
        pdfs = Weights([database["pdf"][prg][ix] for ix in get_idx(prg, ttp, filter=filter)])
    end

    qts = 0.:0.0025:1.
    quantiles = []
    for q in qts 
        push!(quantiles, quantile(vals, pdfs, q))
    end 

    diff_quantiles = diff(quantiles)

    median_quantiles = median(diff_quantiles)
    avg, stdev = mean_and_std(diff_quantiles)
    @printf(" Median Quantile jump: %5.2f  (avg %5.2f - std %5.2f)\n", median_quantiles, avg, stdev)
    
    @printf("%25s  %s\n", "High Quantile Jumps", "dm  at M  (quantile)")

    jumps = 0
    for dq in diff_quantiles 
        if dq > avg + 3*stdev
            ix = find_nearest(diff_quantiles, dq)
            @printf("%25s  %5.2f at %5.2f (%4.1f%%)\n", " ", dq, quantiles[ix], qts[ix]*100)
            jumps += 1
        end
    end 
    if jumps == 0
      @printf("%25s\n", "No Quantile Jumps found")
    end


    @printf("%25s  %s\n", "Top 3 Quantile Jumps", "dm  at M  (quantile)")
    n_jumps = 3
    top_three_jumps = partialsortperm(diff_quantiles, rev=true, 1:n_jumps)
    jumps = 0
    for i in 1:1:n_jumps 
        jump = diff_quantiles[top_three_jumps[i]]
        if jump > avg + 2*stdev
            ix = top_three_jumps[i]
            @printf("%25s  %5.2f at %5.2f (%4.1f%%)\n", " ", jump, quantiles[ix], qts[ix]*100)
            jumps += 1
        end
    end 
    if jumps == 0
      @printf("%25s\n", "No Quantile Jumps found")
    end



end 


function report_model(f, MODELS, logM, q, logP, columns)
    model = MODELS[to_key(logM)*"_"*to_key(logP)*"_"*to_key(q)]
    n = model["mod"]
    print_history(f, n, columns, "2")
    @printf("Star 1 termination:     %s\n", model["end1"])
    @printf("Star 2 termination:     %s\n", model["end2"])
    @printf("Reason for termination: %s\n", model["casus"])
    @printf("Mass transfer case(s):  %s\n", model["case"])
    @printf("\n")
    @printf("SN order : %s - %s\n\n", model["first"], model["second"])
    primary_SN = model["first"]=="1" ? "1" : "2"


    @printf("First SN "); which_SN = primary_SN
    for key in keys(model["SN"*which_SN])
        model["SN"*which_SN][key] > 0 && (@printf("%10s %5s ", "Type", key); break)
    end 
    @printf("solver pathway: %s\n", model["endvals"*which_SN].solver_string)
    @printf("Mf %5.1f Mhe %5.1f Mco %5.1f\n", 
            model["M"*which_SN*"end"], model["M"*which_SN*"_he"], model["M"*which_SN*"_co"])
    @printf("dM_C %.2f Eexp %5.1e  Mni %5.1f Mremg %5.1f\n",
            model["endvals"*which_SN].deltaM_C, 
            model["endvals"*which_SN].E_exp, 
            model["endvals"*which_SN].M_ni, 
            model["endvals"*which_SN].M_remnant_g, )

    secondary_SN = model["first"]=="1" ? "2" : "1"
    @printf("Second SN "); which_SN = secondary_SN
    for key in keys(model["SN"*which_SN])
        model["SN"*which_SN][key] > 0 && (@printf("%10s %5s ", "Type", key); break)
    end 
    @printf("solver pathway: %s\n", model["endvals"*which_SN].solver_string)

    @printf("Mf %5.1f Mhe %5.1f Mco %5.1f\n", 
            model["M"*which_SN*"end"], model["M"*which_SN*"_he"], model["M"*which_SN*"_co"])
    @printf("dM_C %.2f Eexp %5.1e  Mni %5.1f Mremg %5.1f\n",
            model["endvals"*which_SN].deltaM_C, 
            model["endvals"*which_SN].E_exp, 
            model["endvals"*which_SN].M_ni, 
            model["endvals"*which_SN].M_remnant_g, )
end

function print_history2(f, MODELS, logM, q, logP, columns)
    model = MODELS[to_key(logM)*"_"*to_key(logP)*"_"*to_key(q)]
    n = model["mod"]
    print_history(f, n, columns, "2")
end


function report_model_detail(f_MW, MODELS, logM, q, logP, columns, which; details=false, resolve_phase = [10,10,10,10,10], return_arrays = false)
    model = MODELS[to_key(logM)*"_"*to_key(logP)*"_"*to_key(q)]
    n = model["mod"]
    if details
        @printf("Star 1 termination:     %s\n", model["end1"])
        @printf("Star 2 termination:     %s\n", model["end2"])
        @printf("Reason for termination: %s\n", model["casus"])
        @printf("Mass transfer case(s):  %s\n", model["case"])
        @printf("\n")
        @printf("SN order : %s - %s\n\n", model["first"], model["second"])
        primary_SN = model["first"]=="1" ? "1" : "2"


        @printf("First SN "); which_SN = primary_SN
        for key in keys(model["SN"*which_SN])
            model["SN"*which_SN][key] > 0 && (@printf("%10s %5s ", "Type", key); break)
        end 
        @printf("solver pathway: %s\n", model["endvals"*which_SN].solver_string)
        @printf("Mf %5.1f Mhe %5.1f Mco %5.1f\n", 
                model["M"*which_SN*"end"], model["M"*which_SN*"_he"], model["M"*which_SN*"_co"])
        @printf("dM_C %.2f Eexp %5.1e  Mni %5.1f Mremg %5.1f\n",
                model["endvals"*which_SN].deltaM_C, 
                model["endvals"*which_SN].E_exp, 
                model["endvals"*which_SN].M_ni, 
                model["endvals"*which_SN].M_remnant_g, )

        secondary_SN = model["first"]=="1" ? "2" : "1"
        @printf("Second SN "); which_SN = secondary_SN
        for key in keys(model["SN"*which_SN])
            model["SN"*which_SN][key] > 0 && (@printf("%10s %5s ", "Type", key); break)
        end 
        @printf("solver pathway: %s\n", model["endvals"*which_SN].solver_string)

        @printf("Mf %5.1f Mhe %5.1f Mco %5.1f\n", 
                model["M"*which_SN*"end"], model["M"*which_SN*"_he"], model["M"*which_SN*"_co"])
        @printf("dM_C %.2f Eexp %5.1e  Mni %5.1f Mremg %5.1f\n",
                model["endvals"*which_SN].deltaM_C, 
                model["endvals"*which_SN].E_exp, 
                model["endvals"*which_SN].M_ni, 
                model["endvals"*which_SN].M_remnant_g, )
    end

    hf_name = @sprintf("/vol/aibn133/data1/aercolino/GRID_CACHES/hjin/reduced_grid/%5.3f/%5.3f_%5.3f_%s_new.data", logM, q, logP, which)
    hf = CSV.read(hf_name, header=1, DataFrame,delim=' ',ignorerepeated=true, missingstring="NaN") 
    
    age = hf.star_age 
    hf_end = length(age)    
    ZAMS, H_burn, TAMS, end_Hburn, ini_Heburn, He_burn, end_Heburn, ini_Cburn, end_Cburn = evolutionary_checkpoints(hf, hf_end)

    function substager(i,f,d)
        if f > 0
            range_i = []#[i] 
            range_b = d > 0 ? (i:(Int(floor( (f-i)/d ))):f) : []
            range_f = []#d == 0 ? [f] : []
            return union(range_i, range_b, range_f)
        else 
            return []
        end
    end
    substage1 = substager(ZAMS, TAMS, resolve_phase[1])
    substage2 = substager(TAMS, ini_Heburn, resolve_phase[2])
    substage3 = substager(ini_Heburn, end_Heburn, resolve_phase[3])
    substage4 = substager(end_Heburn, ini_Cburn, resolve_phase[4])
    substage5 = substager(end_Heburn, hf_end, resolve_phase[5])
    stages_all = union(substage1,substage2,substage3,substage4, substage5)
    print("\n")
    @printf("Detailed History for model logM-logP-q = %5.3f-%5.3f-%5.3f\n", logM, logP, q)
    print("STAGE           |  Age (Myr) ")
    for col in columns 
        @printf(" %10s", col)
    end 
    print("\n")

    for stage in substage1
        @printf("%16d%6s|", stage, "MS")
        @printf("%12.2f", age[stage]/1e6)
        for col in columns 
            @printf(" %10.2f", hf[!, col][stage])
        end 
        @printf("\n")
    end 
    for stage in substage2
        @printf("%10d%6s|", stage, "pMS")
        @printf("%12.2f", age[stage]/1e6)
        for col in columns 
            @printf(" %10.2f", hf[!, col][stage])
        end 
        @printf("\n")
    end 
    for stage in substage3
        @printf("%10d%6s|", stage, "HeB")
        @printf("%12.2f", age[stage]/1e6)
        for col in columns 
            @printf(" %10.2f", hf[!, col][stage])
        end 
        @printf("\n")
    end 
    for stage in substage4
        @printf("%10d%6s|", stage, "pHeB")
        @printf("%12.2f", age[stage]/1e6)
        for col in columns 
            @printf(" %10.2f", hf[!, col][stage])
        end 
        @printf("\n")
    end 
    for stage in substage5
        @printf("%10d%6s|", stage, "CB")
        @printf("%12.2f", age[stage]/1e6)
        for col in columns 
            @printf(" %10.2f", hf[!, col][stage])
        end 
        @printf("\n")
    end 

    if !return_arrays 
        return nothing 
    end 

    output = Dict() 

    output["star_age"]=hf.star_age[minimum(stages_all):maximum(stages_all)]
    for col in columns 
        output[col]=hf[!, col][minimum(stages_all):maximum(stages_all)]
    end 
    output["label"] = @sprintf("[%1s] logM-logP-q: %5.2f-%5.2f-%5.2f", which, logM, logP, q)
    return output
end


function collect_from_models(;
        val      =(what="orbit", from="post_SN_orbit"),
        condition=(what="orbit", from="post_SN_orbit", operator = ==, threshold=0.0),
        
        MOD=MODELS)
    vals = Float64[]
    pdfs = Float64[]
    
    for model_key in keys(MOD)
        model = MOD[model_key]
        
        
        if from == "post_SN_orbit"
            val = [model["post_SN_orbit"][ix][what] for ix in range(1,length(model["post_SN_orbit"]))]
        elseif from == "endvals"
            val = model["1"].endvals[1]._M_co  # Adjust based on your needs
        elseif from == "SN"
            val = model["1"]["SN"]
        # Add more extraction patterns as needed
        end
        
        push!(vals, val)
        push!(pdfs, model["pdf"]["tot"])
    end
    
    return vals, pdfs
end

function get_from_MODELS(path, model)

    length(path) == 0 && return model
    dict_layer = path[1]
    if dict_layer ∈ ["1", "2", "M", "S", 
                     "1-MT",   "2-MT", "S+NI", "tot" ]
        firstSN = model["1_whichstar"]
        primary   = (firstSN == "1") ? "1" : "2"
        secondary = (firstSN == "1") ? "2" : "1"
        mapper = Dict("M"   =>["1"], 
                      "S"   =>["1"],
                      "1"   =>[primary],
                      "1-MT"=>[primary],
                      "2"   =>[secondary],
                      "2-MT"=>[secondary],
                      "S+NI"=>["1", "2"],
                      "tot" =>["1", "2"])
        postSNruns=length(model["2"])

        for mapped in mapper[dict_layer]
            RLOF_pre1stSN = length(mod["1"].case[1]) > 0 
            RLOF =  if (primary == firstSN) 
                    RLOF_pre1stSN 
                else
                    RLOF_pre1stSN || length(mod["2"].case[s])>0
                end
            okM = ( firstSN == "M" && dict_layer == "M")
            ok1 = ( firstSN != "M" && dict_layer ∈ ["1-MT","2-MT"]         && RLOF)
            ok2 = ( firstSN != "M" && dict_layer ∈ ["1", "2", "S", "S+NI"] && !RLOF)
            ok3 = ( firstSN != "M" && dict_layer ∈ ["1", "2"])
        end

    elseif dict_layer ∈ ["1_whichstar", "casus", 
                         "end_primary", "end_secondary"]
        return model[dict_layer]




        
    elseif dict_layer ∈ ["unboundstate"]
        length(path) == 1 && return model[dict_layer]
        entry = path[2]
        if entry in keys(model[dict_layer]) 
            return model[dict_layer][entry]
        else 
            @printf("ERROR! Wrong path specified\n")
            @printf("The code  attempted to read\nMODELS[\"%s\"][\"%s\"]\n", dict_layer, entry)
            @printf("The available keys are:\n")
            print(keys(model[dict_layer]))
            throw(ErrorException("Terminating: wrong MODEL path (2nd layer)"))
        end
    elseif dict_layer ∈ ["ECI"]
        length(path) == 1 && return model[dict_layer]
        output = []
        for ix in range(1,length(model[dict_layer]))
            if path[2] ∈ ["manual", "ogata"]
                if length(path) == 3
                    push!(output, model[dict_layer][ix][path[2]][path[3]])
                else 
                    return model[dict_layer][ix][path[2]]
                end
            else 
                push!(output,  model[dict_layer][ix][path[2]])
            end
        end
        length(output)>0 && return output
        @printf("ERROR! Wrong path specified\n")
        @printf("The code  attempted to read\nMODELS")
        [@printf("[\"%s\"]", path[i]) for i in range(1,length(path))]
        @printf("\nThe available keys are:\n")
        print(keys(model[dict_layer]))
        throw(ErrorException("Terminating: wrong MODEL path (2nd layer)"))

    elseif dict_layer ∈ ["pdf"]


    elseif dict_layer ∈ ["post_SN_orbit"]

    else
        @printf("ERROR! Wrong path specified\n")
        @printf("The code read MODELS[\"%s\"]\n", dict_layer)
        @printf("The available keys are:\n")
        print(keys(model))
        throw(ErrorException("Terminating: wrong MODEL path (1st layer)"))
    end

end



function seqdata_MODELS(what, names, MODELS; full_output=true, norm_pdf=1)
    @printf("Getting data for %s...\n", what)
    # This function collects data from the models based on the specified 'what' parameter
    # and returns a sequence of data for each model in 'names'. 
    
    full_output && (out_prg= String[] )
    full_output && (out_pdf= Float64[])
    full_output && (out_sn1= Bool[])
    out_what=[]

    for name in names 
        model=MODELS[name]

        firstSN = model["1_whichstar"]

        merger = firstSN == "M"
        single = firstSN == "S"
        primary   = (firstSN in ["1", "2"]) ? (firstSN=="1" ? "1" : "2") : ""
        secondary = (firstSN in ["1", "2"]) ? (firstSN=="1" ? "2" : "1") : ""

        elements=["1"]
        ! (single || merger) && push!(elements, "2")
        ps=Dict(primary=>"1", secondary=>"2")
        n_runs_total = maximum( [length(model[e].SN) for e in elements] )
        weights= (merger || single) ? Dict("ix"=>0,"counter"=>0) : model["unboundstate"]
       
        total_runs = n_runs_total +  max(weights["counter"],1) -1
        pdf = model["pdf"].tot/total_runs
        RLOF_preSN  = length(model["1"].case[1]) > 0 

        for e in elements 
            mod=model[e]
            n_runs = length(model[e].SN)

            full_output && (sub_prg = String[] )
            full_output && (sub_pdf = Float64[])
            full_output && (sub_sn1 = Bool[])
            sub_what= []

            for n in range(1,n_runs)
                
                RLOF = nothing 
                if !(single || merger)
                    RLOF = (e == "1") ? RLOF_preSN : (RLOF_preSN|| (length(model["2"].case[n]) > 0) )
                end

                prg  = nothing
                case = mod.case[n]
                if merger
                    prg = firstSN
                elseif single 
                    prg = "NI-S" 
                else
                    prg = ps[e]
                    prg = RLOF ? "MT-"*prg : "NI-"*prg
                end

                res = get_data_from_model(model, e, n, what)

                full_output && (push!(sub_prg, prg        ))
                full_output && (push!(sub_pdf, pdf        ))
                full_output && (push!(sub_sn1, (e == "1") ))
                push!(sub_what, res)
            end
            if n_runs == 1
                    full_output && (sub_prg  = [sub_prg[1]  for _ in range(1,n_runs_total)])
                    full_output && (sub_pdf  = [sub_pdf[1]  for _ in range(1,n_runs_total)])
                    full_output && (sub_sn1  = [sub_sn1[1]  for _ in range(1,n_runs_total)])
                    for n in 2:n_runs_total 
                        res = get_data_from_model(model, e, n, what)
                        push!(sub_what, res)
                    end
            end

            (full_output && (weights["ix"] > 0)) && (sub_pdf[weights["ix"]] *= max(weights["counter"]))

            
            for i in range(1,length(sub_what))
                full_output && (push!(out_prg, sub_prg[i]))
                full_output && (push!(out_pdf, sub_pdf[i]))
                full_output && (push!(out_sn1, sub_sn1[i]))
                push!(out_what, sub_what[i])
            end
        end
    end
    full_output && (out_pdf = out_pdf ./ norm_pdf)
    full_output && return out_prg, out_pdf, out_sn1, out_what
    return out_what
end

function get_data_from_model(model, which, n_from_call, what::String)
    mod=model[which]
    # throw(ErrorException)
    data_from_eci = occursin("ECI?", what) || occursin("Rmax", what) || what == "tau_ECI" || what ∈ ["Eheat", "mh", "Omega_Eff", "E_B"] 
    data_from_postSN =  what ∈ ["a", "a_peri", "P", "ϵ", "e", 
                                "orbit", "will_one_periastron_occur",
                                "when_will_periastron_occur" ]
    n=n_from_call
    if which == "1"
        n = (data_from_eci || data_from_postSN || what == "1st_kick") ? n_from_call : 1
    end

    what=="SN" && return mod.SN[n]
    if what == "1st_kick"
        try
            return model["1st_kick"][n]
        catch 
            print(keys(model))
        end
    end
    what ∈ ["Mmax", "dM_C", "t_SN"] && return getproperty(mod, Symbol(what))[n] 
    what ∈ ["M_end", "M_he", "M_co", "Xc",
            "deltaM_C", "M_ni", "E_exp", 
            "M_remnant_g", "M_remnant_b",
            "v_kick", "solver_string"] && return getproperty(mod.endvals[n], Symbol(what)) 
    
    data_from_postSN && return getproperty(model["post_SN_orbit"][n], Symbol(what))
    
    if what ∈ ["pre_v1", "pre_a", "pre_v2"]
        try
            return getproperty(model["pre_SN_orbit"], Symbol(what[5:end]))
        catch 
            print(keys(model))
        end
    end

    if what == "kick_angle"
        # print(keys(model))
        x, y, z = model["1st_kick_whole"][n]
        r = sqrt(x^2+y^2+z^2)
        φ = asin(z/r)
        θ = atan(x/(r*cos(φ)), y/(r*cos(φ)))
        return (θ/π,φ/π)
    end

    if what == "1stCO"
        outcome = model["1"].SN[1] 
        explosion = outcome in [:IIP, :IIn, :SN87A, :IIb, :Ibc, :Ibn]
        explosion && return :NS 
        outcome == :BH && return :BH 
        outcome == :WD && return :WD 
        throw(ErrorException("Failed to assign a 1st compact-object: outcome = $outcome; n=$n"))
    end
    if what == "2ndCO"
        outcome = model["2"].SN[n] 
        explosion = outcome in [:IIP, :IIn, :SN87A, :IIb, :Ibc, :Ibn]
        explosion && return :NS 
        outcome == :BH && return :BH 
        outcome == :WD && return :WD 
        throw(ErrorException("Failed to assign a 1st compact-object: outcome = $outcome; n=$n"))
    end

    what == "mod" && return model["mod"]


    eci=model["ECI"][n]

    occursin("ECI?", what) && return eci[what[6:end]]["ECI?"] 
    what == "tau_ECI" && return eci["ogata"][what] 
    occursin("Rmax", what) && return eci[what[6:end]]["Rmax"] 
    what ∈ ["Eheat", "mh", "Omega_Eff", "E_B"] && return eci[what]  
    
    throw(ErrorException("Undefined \"what\" in get_data_from_model. Got $what.\n"))
    
end
