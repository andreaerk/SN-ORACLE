using TOML 
show_uuid(package_name::AbstractString) = TOML.parsefile(Base.active_project())["deps"][package_name]
