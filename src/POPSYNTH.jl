module POPSYNTH 

    using CSV, DataFrames, DataStructures
    using Printf, Debugger
    using Meshes
    using Interpolations 
    using LinearAlgebra
    using QuadGK
    using Statistics 
    using StatsBase
    using NearestNeighbors
    using StatsBase 
    using HDF5 
    using Dierckx
    using Distributions 

    const workdir = "/vol/aibn133/data1/aercolino/SOFTWARE/SN_from_grids/" 
    const out_dir = workdir * "output/"
    const data_dir = workdir * "Data/" 
    const pythondir = workdir * "PythonScripts/"
    const SG_dir = data_dir * "Models_SG/"
    const BG_dir = data_dir * "Models_BG/"
    const CC_dir = data_dir * "Models_CC/"
    const S24_dir = data_dir * "Schneider2024/"
    const PS20_dir = data_dir * "PS2020/"

    include("defs.jl")
    include("kicks.jl")
    include("CaseBB.jl")
    include("SN_support.jl")
    include("SN_reader_MW.jl")
    include("SN_data_read.jl")

    include("ECI.jl")
    using .ECI

    export predict_SN_from_each_model, do_SN_popsynth, reorganize_SN_data, get_data_sne
    export G, Msun, Rsun, yr, kB, km, kms, day, Lsun
    
end #module POPSYNTH