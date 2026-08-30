function apply_kick_and_recenter(m_rem, v1, m2, v2; v_kick=[0.,0.,0.])

    #assume v1∥-v2 -> coplanar 
    v_1 = [0, +v1, 0]
    v_2 = [0, -v2, 0] 

    v_1_postSN = v_1 + v_kick
    v_2_postSN = v_2 
    v_com_postSN = (v_2_postSN * m2 + v_1_postSN * m_rem) / (m2 + m_rem)

    v_1_new = v_1_postSN - v_com_postSN
    v_2_new = v_2_postSN - v_com_postSN

    return (v_1_new, v_2_new)
end

function primitive_freefall_t(r, μ,  ϵ,  infall)
    sign = infall ? -1 : +1

    c = sign * 1/sqrt(2)
    t1 = r/ϵ * sqrt(ϵ+ μ/r) 
    t2 = μ/sqrt(abs(ϵ)^3) * atan(sqrt(ϵ+ μ/r) / sqrt(abs(ϵ)))
    return c * (t1-t2)
end


function derive_orbital_solution(M1, R1, V1, M2, R2, V2, impact_distance; G = G)
    t = 0


    μ = G * (M1 + M2) #gravitational parameter

    #relative distance vector
    r_vec = R2 - R1 
    r = norm(r_vec) 
    r_hat = r_vec / r
    #relative velocity vector
    v_vec = V2 - V1 
    v = norm(v_vec) 
    v_hat = v_vec / v
    
    # current configuration is getting closer or further away from 0?
    local_rdot = r_hat ⋅ v_vec
    approaxing_periaxis = local_rdot < 0

    #specific energy
    ϵ = 1/2 * v^2 - μ/r 
    
    #is orbit bound?
    orbit = (ϵ < 0) ? (:bound) : (:unbound)

    #specific angular momentum
    h_vec = r_vec × v_vec
    h = norm(h_vec)
    h_hat = h_vec / h 

    #if r∥v -> direct collision or escape
    if abs(h) < 1e-12 
        if r_vec ⋅ v_vec < 0 || ϵ < 0 #direct collision
            t_infall = primitive_freefall_t(impact_distance, μ,  ϵ,  true)-primitive_freefall_t(r, μ,  ϵ,  true)
            
            if r_vec ⋅ v_vec > 0
                r_max = - μ/ϵ
                t_infall += 2*(primitive_freefall_t(r_max, μ,  ϵ,  false)-primitive_freefall_t(r, μ,  ϵ,  false))
            end
            if t_inflall < 0 
                @printf("WARNING: negative time to infall! r=%f, v=%f, μ=%f, ϵ=%f\n", r, v, μ, ϵ)
                throw(ErrorException)
            end
            return (a=NaN, a_peri=0,
            P=NaN, e=NaN, ϵ=ϵ, 
            orbit = :direct_collision, 
            will_one_periastron_occur=true, 
            when_will_periastron_occur = t_infall,
            )
        else #escape
            return (a=NaN, a_peri=r,
            P=NaN, e=NaN, ϵ=ϵ, 
            orbit = :unbound, 
            will_one_periastron_occur=false, 
            when_will_periastron_occur = 0)
        end
    end

    #eccentricity vector
    e_vec = v_vec × h_vec / μ - r_hat   
    e = norm(e_vec)
    e_hat = e_vec / e
    
    #p-vector
    p_vec = h_vec × e_vec
    p = norm(p_vec) 
    p_hat = p_vec / p 

    #periastron calculation
    a_periastron = h^2 / (μ * (1+e))
    #ν = true anomaly
    cos_ν = e_hat ⋅ r_hat
    sin_ν = p_hat ⋅ r_hat

    ν = atan( sin_ν , cos_ν) #defined between [-π, +π]

    #will one periastron occur? if orbit is bound (ϵ<0 || e<1) or
    #                           if periapsis is yet to occur (sin(ν)<0)
    will_one_periastron_occur = (e < 1) || approaxing_periaxis 
    
    #ν ∈ [0, 2π]
    # ν = atan( sin_ν / cos_ν) + π * (sin_ν < 0 ?  (cos_ν >= 0 ? +2 : +1) : (cos_ν >= 0 ? 0 : 1))
    
    if will_one_periastron_occur && orbit in [:bound,:unbound]
        M_eff = nothing
        if abs(1-e) < 1e-12  #parabolic orbit
            # Barker formula
            D = sign(ν) * tan(ν/2) 
            p = h^2/μ
            t = 0.5 * sqrt(p^3/(2μ)) * (D + D^3/3)  # signed: >0 future, <0 past
            if t < 0 
                @printf("WARNING: negative time to periastron (parabola)\n")
                @printf("     r = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", r, r_vec[1], r_vec[2], r_vec[3])
                @printf("     v = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", v, v_vec[1], v_vec[2], v_vec[3])
                @printf("     h = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", h, h_vec[1], h_vec[2], h_vec[3])
                @printf("     e = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", e, e_vec[1], e_vec[2], e_vec[3])
                @printf("     μ = %7.2e\n", μ)
                @printf("     ϵ = %7.2e\n",ϵ)
                @printf("     ν = %7.2f, cos = %7.2f, sin = %7.2f\n", ν, cos_ν, sin_ν)
                @printf("     D = %7.2f\n",    D)
                @printf("     p = %7.2f\n",    p)
                @printf("     t = %7.2e\n",    t)
                throw(ErrorException)
            end

        elseif e > 1 #hyperbolic orbit
            if abs(1 + e*cos_ν) < 1e-10
                @printf("warning - hyperbolic motion close to parabolic!")
            end
            coshF = (e + cos_ν) / (1 + e*cos_ν)
            coshF = max(coshF, 1.0)             # numeric safety
            F = acosh(coshF)          # signed hyperbolic anomaly
            M = e*sinh(F) - F                   # mean hyperbolic anomaly
            # time scale for hyperbola: sqrt((-a)^3/μ)
            a = -μ/(2*ϵ)                         # negative for hyperbola
            t = sqrt((-a)^3/μ) * M     # t<0 -> past, t>0 -> future
            if t < 0 
                @printf("WARNING: negative time to periastron (hyperbola)\n")
                @printf("     r = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", r, r_vec[1], r_vec[2], r_vec[3])
                @printf("     v = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", v, v_vec[1], v_vec[2], v_vec[3])
                @printf("     h = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", h, h_vec[1], h_vec[2], h_vec[3])
                @printf("     e = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", e, e_vec[1], e_vec[2], e_vec[3])
                @printf("     μ = %7.2e\n", μ)
                @printf("     ϵ = %7.2e\n",ϵ)
                @printf("     ν = %7.2f, cos = %7.2f, sin = %7.2f\n", ν, cos_ν, sin_ν)
                @printf("cosh F = %7.2f\n",coshF)
                @printf("     F = %7.2f\n",    F)
                @printf("     M = %7.2f\n",    M)
                @printf("     a = %7.2e\n",    a)
                @printf("     t = %7.2e\n",    t)
                throw(ErrorException)
            end

        else #elliptic orbit
            if abs(1 + e*cos_ν) < 1e-16
                @printf("warning - high eccentricity!")
            end
            cosE = (e + cos_ν) / (1 + e*cos_ν)
            sinE = sqrt(max(0.0, 1 - e^2)) * sin_ν / (1 + e*cos_ν)
            E = atan(sinE, cosE)
            M = E - e*sin(E)
            a = -μ/(2*ϵ)
            t = M / sqrt(μ/a^3)         
            # Ensure t is time to next periastron passage (always positive)
            if t < 0
                a = - μ/(2*ϵ)
                P = 2π*sqrt(a^3/μ)
                t += P
                if t < 0 
                    @printf("WARNING: negative time to periastron (elliptic)\n")
                    @printf("    r = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", r, r_vec[1], r_vec[2], r_vec[3])
                    @printf("    v = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", v, v_vec[1], v_vec[2], v_vec[3])
                    @printf("    h = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", h, h_vec[1], h_vec[2], h_vec[3])
                    @printf("    e = %7.2e  ⃗[%7.2e %7.2e %7.2e]\n", e, e_vec[1], e_vec[2], e_vec[3])
                    @printf("    μ = %7.2e\n", μ)
                    @printf("    ϵ = %7.2e\n",ϵ)
                    @printf("    ν = %7.2f, cos = %7.2f, sin = %7.2f\n", ν, cos_ν, sin_ν)
                    @printf("cos E = %7.2f\n",cosE)
                    @printf("    E = %7.2f\n",   E)
                    @printf("    M = %7.2f\n",   M)
                    @printf("    a = %7.2e\n",   a)
                    @printf("    t = %7.2e\n",   t)
                    throw(ErrorException)
                end

            end
        end
    end
    #semimajor axis
    a = (ϵ<0) ? - μ/(2*ϵ) : NaN
    #orbital period
    P = (ϵ<0) ? 2π*sqrt(a^3/μ) : NaN
    
    a_periastron <= impact_distance && (orbit = :direct_collision)

    return (a=a, a_peri=a_periastron, P=P, e=e, ϵ=ϵ, 
            orbit = orbit, 
            will_one_periastron_occur=will_one_periastron_occur, 
            when_will_periastron_occur = t)
end




function test_orbital_function()
    
    cases = [
        (
            "Circular orbit, radius 1",
            1.0, [-0.5, 0.0, 0.0], [0.0, -sqrt(2)/2, 0.0],
            1.0, [+0.5, 0.0, 0.0], [0.0, +sqrt(2)/2, 0.0], 1
        ),       
        (
            "Elliptical e=0.5, a=2",
            1.0, [0.0, 0.0, 0.0], [0.0, 0.5, 0.0],
            1.0, [1.0, 0.0, 0.0], [0.0, 1.5, 0.0], 1
        ),
        (
            "Hyperbolic e=1.5",
            1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
            1.0, [5.0, 0.0, 0.0], [-0.5, 1.2, 0.0], 1
        ),
        (
            "Direct collision (1)",
            1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
            1.0, [5.0, 0.0, 0.0], [-.7, 0.0, 0.0], 1
        ),
        (
            "Direct collision (2)",
            1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
            1.0, [5.0, 0.0, 0.0], [+.7, 0.0, 0.0], 1
        ),
        (
            "Direct escape",
            1.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
            1.0, [5.0, 0.0, 0.0], [1.5, 0.0, 0.0], 1
        ),
        (
            "Sun-Earth",
            1*Msun, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
            3.00274e-6*Msun, [0.983*AU, 0.0, 0.0], [0, 30.29*km/s, 0.0], G
        ),
        (
            "Earth-Moon",
            3.00274e-6*Msun, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
            3.6943e-8*Msun, [405500*km, 0.0, 0.0], [0, 0.970*km/s, 0.0], G
        ),
        # Elliptic (a=2, e=0.5), μ = G*(1+1) = 2
        (
        "Elliptic (a=2,e=0.5) — before periapsis",
        1.0, [-0.45325421887794337,  0.2616864452805141, 0.0], [-0.28867513459481287, -0.7886751345948129, 0.0],
        1.0, [ 0.45325421887794337, -0.2616864452805141, 0.0], [ 0.28867513459481287,  0.7886751345948129, 0.0], 1.0
        ),
        (
        "Elliptic (a=2,e=0.5) — at periapsis",
        1.0, [-0.5, 0.0, 0.0], [ 0.0, -0.8660254037844386, 0.0],
        1.0, [ 0.5, 0.0, 0.0], [ 0.0,  0.8660254037844386, 0.0], 1.0
        ),
        (
        "Elliptic (a=2,e=0.5) — after periapsis",
        1.0, [-0.45325421887794337, -0.2616864452805141, 0.0], [ 0.28867513459481287, -0.7886751345948129, 0.0],
        1.0, [ 0.45325421887794337,  0.2616864452805141, 0.0], [-0.28867513459481287,  0.7886751345948129, 0.0], 1.0
        ),

        # Parabolic (p=2, e=1), μ = 2
        (
        "Parabolic (p=2,e=1) — before periapsis",
        1.0, [-0.4641016151377546,  0.26794919243112264, 0.0], [-0.25, -0.9330127018922194, 0.0],
        1.0, [ 0.4641016151377546, -0.26794919243112264, 0.0], [ 0.25,  0.9330127018922194, 0.0], 1.0
        ),
        (
        "Parabolic (p=2,e=1) — at periapsis",
        1.0, [-0.5, 0.0, 0.0], [ 0.0, -1.0, 0.0],
        1.0, [ 0.5, 0.0, 0.0], [ 0.0,  1.0, 0.0], 1.0
        ),
        (
        "Parabolic (p=2,e=1) — after periapsis",
        1.0, [-0.4641016151377546, -0.26794919243112264, 0.0], [ 0.25, -0.9330127018922194, 0.0],
        1.0, [ 0.4641016151377546,  0.26794919243112264, 0.0], [-0.25,  0.9330127018922194, 0.0], 1.0
        ),

        # Hyperbolic (e=1.5, p=3), μ = 2
        (
        "Hyperbolic (e=1.5,p=3) — before periapsis",
        1.0, [-0.5650354826521339,  0.32622338801089956, 0.0], [-0.20412414523193156, -0.9659258262890683, 0.0],
        1.0, [ 0.5650354826521339, -0.32622338801089956, 0.0], [ 0.20412414523193156,  0.9659258262890683, 0.0], 1.0
        ),
        (
        "Hyperbolic (e=1.5,p=3) — at periapsis",
        1.0, [-0.6, 0.0, 0.0], [ 0.0, -1.0206207261596576, 0.0],
        1.0, [ 0.6, 0.0, 0.0], [ 0.0,  1.0206207261596576, 0.0], 1.0
        ),
        (
        "Hyperbolic (e=1.5,p=3) — after periapsis",
        1.0, [-0.5650354826521339, -0.32622338801089956, 0.0], [ 0.20412414523193156, -0.9659258262890683, 0.0],
        1.0, [ 0.5650354826521339,  0.32622338801089956, 0.0], [-0.20412414523193156,  0.9659258262890683, 0.0], 1.0
        )

    ]

    for (name, m1, r1, v1, m2, r2, v2, G) in cases
        println("=== $name ===")
        res = derive_orbital_solution(m1, r1, v1, m2, r2, v2, G= G)
        @printf("a = %9.2e, a_peri = %9.2e, P = %11.4e, e=%5.3f, ε=%9.2e, orbit? %20s, periastron_passage? %6s if so, when? %9.2e\n", res.a, res.a_peri, res.P, res.e, res.ϵ, res.orbit, res.will_one_periastron_occur, res.when_will_periastron_occur )
    end
end

function test_orbital_function_new()
    G = 1.0

    # 1) Circular, equal masses, separation = 1
    m1, m2 = 1.0, 1.0
    μ = G*(m1+m2)
    r1 = [-0.5, 0.0, 0.0]; r2 = [0.5, 0.0, 0.0]   # separation 1
    v_rel = sqrt(μ / 1.0)                         # relative circular speed
    v1 = [0.0, -v_rel/2, 0.0]; v2 = [0.0, v_rel/2, 0.0]

    # 2) Elliptical example: choose a=2, e=0.5 and place at periastron (r=a(1-e)=1)
    m1b, m2b = 1.0, 1.0
    μb = G*(m1b+m2b)
    a = 2.0; e = 0.5
    rp = a*(1-e) # = 1
    r1b = [-rp/2, 0.0, 0.0]; r2b = [rp/2, 0.0, 0.0]
    h = sqrt(μb * a * (1 - e^2))
    v_rel_b = h / rp
    v1b = [0.0, -v_rel_b/2, 0.0]; v2b = [0.0, v_rel_b/2, 0.0]

    # 3) Hyperbolic example (your earlier example)
    m1c, m2c = 1.0, 1.0
    rc1 = [0.0, 0.0, 0.0]; vc1 = [0.0, 0.0, 0.0]
    rc2 = [5.0, 0.0, 0.0]; vc2 = [0., 1.2, 0.0]

    println("=== Circular (corrected) ===")
    println(derive_orbital_solution(m1, r1, v1, m2, r2, v2; G=G))
    println("=== Elliptic (a=2, e=0.5 at peri) ===")
    println(derive_orbital_solution(m1b, r1b, v1b, m2b, r2b, v2b; G=G))
    println("=== Hyperbolic example ===")
    println(derive_orbital_solution(m1c, rc1, vc1, m2c, rc2, vc2; G=G))
end



avg_v_kick_Hobbs = 265 
v_kicks_Hobbs() = avg_v_kick_Hobbs * rand(Chi(3))
v_kicks_DM25() = rand(LogNormal(5.6, 0.68))
function v_kicks_COMBINE(stage)
    avg_v_kick=nothing
    if stage == "ZAMS"
        avg_v_kick= avg_v_kick_Hobbs
    elseif stage == "stripped"
        avg_v_kick = (rand() >= 0.8) ? 200 : 120 
    elseif stage == "CaseBB"
        avg_v_kick = (rand() >= 0.8) ? 200 : 60
    elseif stage == "ultrastripped" 
        avg_v_kick = (rand() >= 0.8) ? 200 : 30 
    else
        throw(ErrorException("Unknown COMBINE kick end-stage: $stage"))
    end
    return avg_v_kick * rand(Chi(3)) 
end

function v_kicks_COMBINE2(stage)
    avg_v_kick=nothing
    if stage == "ZAMS"
        return v_kicks_DM25()
    elseif stage == "stripped"
        avg_v_kick = (rand() >= 0.8) ? 200  : 120
        return avg_v_kick * rand(Chi(3)) 
    elseif stage == "CaseBB"
        avg_v_kick = (rand() >= 0.8) ? 200 : 60
        return avg_v_kick * rand(Chi(3)) 
    elseif stage == "ultrastripped" 
        avg_v_kick = (rand() >= 0.8) ? 200 : 30 
        return avg_v_kick * rand(Chi(3)) 
    else
        throw(ErrorException("Unknown COMBINE (2) kick end-stage: $stage"))
    end
end

function v_kicks_Valli25(stage)
    avg_v_kick=nothing
    if stage == "ZAMS"
        avg_v_kick = 300
        return avg_v_kick * rand(Chi(3)) 
    elseif stage == "stripped"
        return rand(Normal(100, 11))
    elseif stage == "CaseBB" 
        avg_v_kick = 5
        return avg_v_kick * rand(Chi(3)) 
    else
        throw(ErrorException("Unknown Valli kick end-stage: $stage"))
    end
end


function v_kicks_MM20(parameter; v_ns=400) #parameter = (Mco-Mns)/Mn 
    return rand(Normal(v_ns*parameter))
end

function get_kick_parameter(KICKS, endvals, case)
    if occursin("MM", KICKS)
        return (endvals.M_co-endvals.M_remnant_b)/endvals.M_remnant_b 
    elseif occursin("Hobbs", KICKS) || occursin("DM25", KICKS)
        return nothing
    elseif occursin("COMBINE", KICKS) || occursin("Valli25", KICKS)
        return case
    elseif KICKS == "Inf" || KICKS == "None" || KICKS[1:2]=="c_"
        return nothing
    else 
        throw(ErrorException("Unknown kick type: $KICKS"))
    end
end



function draw_φ(KICKS, parameter)
    sinφ = 0
    
    if ! occursin("Valli25", KICKS) || parameter != "stripped"
        sinφ = 2*rand() - 1
    elseif occursin("Valli25", KICKS) && parameter == "stripped"
        sinφ =  rand() * ( sin(deg2rad(90)) - sin(deg2rad(85))) + sin(deg2rad(85))
        rand(Bool) && (sinφ *= -1)
    end
    
    return asin(sinφ)
end

function get_kicks(KICKS, parameter, skip_n) #kicks to be given in kms
    num = nothing
    v_kick() = if occursin("MM_", KICKS) 
                v_ns =  if occursin(" _v", KICKS)
                            vk = parse(Int, match(r"_v(\d+)", KICKS).captures[1])
                        else    
                            vk = 400
                        end
                v_kicks_MM20(parameter; v_ns=vk)
              elseif occursin("Hobbs_", KICKS)
                v_kicks_Hobbs()
              elseif occursin("COMBINE_", KICKS)
                v_kicks_COMBINE(parameter)
              elseif occursin("COMBINE2_", KICKS)
                v_kicks_COMBINE2(parameter)
                # throw(ErrorException)
              elseif occursin("DM25_", KICKS)
                v_kicks_DM25()
                elseif occursin("Valli25_", KICKS)
                v_kicks_Valli25(parameter)
              elseif KICKS == "Inf"
                +Inf
              elseif KICKS[1:2] == "c_"
                c_light/kms
              elseif KICKS == "None"
                0
              else 
                throw(ErrorException("Unknown kick type: $KICKS"))
              end
    if occursin("_n", KICKS)
        num = parse(Int, match(r"_n(\d+)", KICKS).captures[1])
    elseif KICKS == "Inf" ||  KICKS == "None"
        num = 1
    else 
        throw(ErrorException("You specified a mode that must be iterated but you don't provide the number of samples. Got $KICKS"))
    end

    postkick_states = []

    for n in range(1, !skip_n ? num : 1)
        # throw(ErrorException)
        θ = rand() * 2π
        φ = draw_φ(KICKS, parameter)
        norm_v = v_kick()
        v = norm_v .* [ cos(θ)*cos(φ), 
                        sin(θ)*cos(φ), 
                               sin(φ) ]
        abs(norm(v) - norm_v)/norm_v > 1e-10 && throw(ErrorException("OH OH!"))
        push!(postkick_states, v)
    end
    

    return postkick_states 
end