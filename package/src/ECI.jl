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




#ECI - Hirai18, Ogata21, Hirai23

module ECI

    G = 6.6743e-8 #g-1 cm3 s-1
    Msun = 1.989e33 #g
    Rsun = 6.955e10 #cm
    c = 2.99792458e10 #cm/s
    yr = 3.154e+7 #s


    Omega_eff(R2, a) = 1/2 * ( 1 - sqrt( 1- (R2/a)^2 ) )
    Eheat(Eexp, OmegaEff; p = 0.10) = p*Eexp*OmegaEff
    Eheat(Eexp, R2,    a; p = 0.10) = p*Eexp*Omega_eff(R2,a)
    eval_mh(Mej, OmegaEff) = 1/2 *  Mej * OmegaEff 
    eval_mh(Mej, R2, a)    = 1/2 *  Mej * Omega_eff(R2, a)
    heating_rate(m, Eheat, mh, M2) = Eheat/mh * min(1, mh/m) / (1 + log(M2/mh))
    heating_rate_mh(Eheat, mh, M2)   = Eheat / (1 + log(M2/mh))
    heating_rate_rest(Eheat, mh, M2) = Eheat * log(M2/mh)/ (1 + log(M2/mh))

    kappa_fit(M2) = 1.24 * (1 - 0.02 * M2)  #M2 in M⊙, result in cm²g⁻¹
    Lmax_analytic(M2) = 4 * pi * G * (M2*Msun) * c / kappa_fit(M2)
    logRmax(Eheat) = -1/14 * (log10(Eheat/0.08e51))^2 + 3.1
    tau_inflation(Eheat, M2; alpha = 0.181) = alpha * Eheat / Lmax_analytic(M2)
    
    function manual_Rmax(Eheat, mh, R2, M2)
        Eb = 0.5 * G * mh*Msun * M2*Msun / (R2 * Rsun) 
        (Eheat >= Eb) && return +Inf 
        return R2 / (1 - Eheat/Eb)  
    end 


    function estimate_ECI(;  Mej::AbstractFloat  = nothing, #Msun
                            Eexp::Float64        = nothing, #erg
                            a::AbstractFloat     = nothing, #Rsun
                            M2::AbstractFloat    = nothing, #Msun 
                            R2::AbstractFloat    = nothing, #Rsun
                            RL2::AbstractFloat   = nothing, #Rsun
                            p::Float64     = 0.08,
                            alpha::Float64 = 0.181)

        Omega_Eff = Omega_eff(R2,a)
        mh = eval_mh(Mej, Omega_Eff)
        Eh = Eheat(Eexp, Omega_Eff, p = p )

        Rmax = max(10. .^ logRmax(Eh), R2)
        tau_ECI = tau_inflation(Eh, M2; alpha=alpha)
        stats_ECI_Ogata = Dict("ECI?"     => (Rmax >= RL2), 
                            "tau_ECI" => tau_ECI,
                            "Rmax"    => Rmax)

        Rmax_manual = max(manual_Rmax(Eh, mh,  R2, M2), R2)  
        stats_ECI_manual = Dict("ECI?" => (Rmax_manual >= RL2), 
                                "Rmax"=> Rmax_manual)

        return Dict("ogata"=> stats_ECI_Ogata, "manual"=>stats_ECI_manual, 
                    "Eheat" => Eh, "mh"=> mh, "Omega_Eff"=>Omega_Eff, "E_B"=> (0.5  * G * M2 * mh *Msun*Msun/ (R2*Rsun)))
    end 


end 

