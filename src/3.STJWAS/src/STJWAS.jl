module ST

using Distributions
using DataFrames
using ..PedModule
using ProgressMeter

include("solver.jl")
include("types.jl")
include("functions.jl")
include("1.markers/markers.jl")
include("2.MCMC/MCMC.jl")
include("variance.jl")
include("output.jl")
include("interface.jl")

export build_model,set_covariate,set_random,add_markers,get_pedigree
export showMME
export outputMCMCsamples,solve,runMCMC

end
