library(JuliaCall)
setwd("C:/Users/jgpen/Documents/Finite-Populations-/Julia-1.2.0/bin")
julia <- julia_setup(JULIA_HOME = "C:/Users/jgpen/Documents/Finite-Populations-/Julia-1.2.0/bin")


julia_install_package_if_needed("Compat")
julia_install_package_if_needed("Statistics")
julia_install_package_if_needed("Random")
julia_install_package_if_needed("Distributions")
julia_install_package_if_needed("DelimitedFiles")
julia_install_package_if_needed("DataFrames")

julia_source("C:/Users/jgpen/Documents/Finite-Populations-/scripts/helpers.jl")
julia_source("C:/Users/jgpen/Documents/Finite-Populations-/scripts/Normal_DMP.jl")
julia_command("model = STotales"); model <- julia_eval("model")


