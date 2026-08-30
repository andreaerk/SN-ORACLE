using Pkg

# Activate POPSYNTH's environment
Pkg.activate(@__DIR__)
Pkg.resolve()

# List of all required packages (external packages and StdLibs)
deps = [
    # External packages
    "CSV", "DataFrames", "DataStructures", "Debugger",
    "Meshes", "Interpolations", "QuadGK",
    "StatsBase", "NearestNeighbors", "HDF5",
    "Dierckx", "Distributions",

    # StdLib packages (need to be in deps explicitly)
    "Printf", "LinearAlgebra", "Statistics"
]

# Add each dependency (no harm if it’s already present)
for dep in deps
    try
        println("Adding dependency: $dep")
        Pkg.add(dep)
    catch e
        @warn "Failed to add $dep" exception=(e, catch_backtrace())
    end
end

# Ensure all dependencies are installed and resolved
Pkg.instantiate()
Pkg.resolve()
println("Dependencies populated and environment instantiated successfully.")
