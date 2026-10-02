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


module SN_ORACLE

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

    script_dir = dirname(@__FILE__)
    const workdir = script_dir * "/../../"
    #const out_dir = workdir * "output/"
    const data_dir = workdir * "Data/" 
    #const pythondir = workdir * "PythonScripts/"
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
    export G, Msun, Rsun, Lsun, yr, kB, km, kms, day

end #module SN_ORACLE
