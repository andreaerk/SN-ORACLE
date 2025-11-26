#!/users/aercolino/.julia/juliaup/julia-1.11.3+0.x64.linux.gnu/bin/julia

include("SN_reader_MW.jl")

#For the Ertl Method, add '-' followed by any of these substrings
# they represent different calibrations for the Ertl method 
#        "S19.8_old" "S19.8_new"
#        "W15_old"   "W15_new"   
#        "W18_old"   "W18_new"   
#        "W20_old"   "W20_new"   
#        "N20_old"   "N20_new"   
#For the Muller Method, add '-' followed by either "M16", "S21" or "AD23" 
# these represent different parameter calibrations 



MERGER_MASSLOSS  = "X"  #"FIX"+fraction or "Energy"
MERGER_EOL       = "CORE"  #"TOTAL" or "CORE" or "ANALYTIC"
fraction_binary = 0.75
KICKS = Inf 

global dummy = Dict()
MODELS = nothing 
# for EXP_CRIT in ["X", "MM-M16", "MM-S21", "MM-AD23", 
#                  "Ertl-S19.8_new", "Ertl-W20_new", "Ertl-W18_new", "Ertl-W15_new", "Ertl-N20_new", 
#                  "PS20", "comp0.20", "comp0.25", "comp0.30", "comp0.35"]
for EXP_CRIT in ["PS20"]
    for MER_CRIT in ["X", "ERK", "PA_IV", "PAULI", "PABLO", "ALL"]

        @printf("importing data... assuming Mer.Crit. %10s   and    Exp.Crit. %10s\n", MER_CRIT, EXP_CRIT)
        

        MODELS = predict_SN_from_each_model(
            MER_CRIT,MERGER_EOL,MERGER_MASSLOSS, EXP_CRIT, KICKS
        )
        ss = [ [alpha_m, 0, 0],
            [alpha_m, alpha_p, alpha_q],
            [alpha_m, alpha_p+stdev_p, alpha_q],
            [alpha_m, alpha_p-stdev_p, alpha_q],
            [alpha_m, alpha_p, alpha_q+stdev_q],
            [alpha_m, alpha_p, alpha_q-stdev_q],
            ]  
        ss = [ [alpha_m, 0, 0] ]
        for s in ss
            sm = s[1]
            sp = s[2]
            sq = s[3]
            for fB in [0.75]#[0, 0.25, 0.50, 0.75, 1.00 ] #
                  _, _, _ = do_SN_popsynth(MER_CRIT, MERGER_EOL, MERGER_MASSLOSS, EXP_CRIT, sm, sp, sq, :custom, fB, MODELS=MODELS)
            end
        end
    end
end

 using Profile
using ProfileView

if true 
    EXP_CRIT="MM"
    MER_CRIT="X"
    MODELS = predict_SN_from_each_model(
        MER_CRIT,MERGER_EOL,MERGER_MASSLOSS, EXP_CRIT, KICKS
    )
    s =   [alpha_m, 0, 0] 
    sm = s[1]
    sp = s[2]
    sq = s[3]
    fB = 0.75
    Profile.clear()

    do_SN_popsynth(MER_CRIT, MERGER_EOL, MERGER_MASSLOSS, EXP_CRIT, sm, sp, sq, :custom, fB, MODELS=MODELS)
    ProfileView.@profview (for i=1:15; _,_,_=  do_SN_popsynth(MER_CRIT, MERGER_EOL, MERGER_MASSLOSS, EXP_CRIT, sm, sp, sq, :custom, fB, MODELS=MODELS); end)

end 