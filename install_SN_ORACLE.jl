using Pkg
using Printf

#developing package
println("--------------------\nDeveloping SN_ORACLE\n--------------------\n")
@printf("%10s: %20s\n%10s: %20s\n", "Author", "Andrea Ercolino", "Contact", "aercolino.astro@gmail.com")
@printf("Latest versions of the package available on\nhttps://github.com/andreaerk/SN-ORACLE/\n\n")
Pkg.develop(path="package")

# Ensure all dependencies are installed and resolved
println("--------------\nInstantiating SN_ORACLE\n--------------\n")
Pkg.instantiate()
println("--------------\nResolving packages for SN_ORACLE\n--------------\n")
Pkg.resolve()

println("\n\nDependencies populated and environment instantiated successfully.\n\n")
