module MT

using Distributions
using DataFrames
using ..PedModule
using ProgressMeter

include("types.jl")
include("solver.jl")
include("functions.jl")
include("residual.jl")
include("MCMC/MCMC.jl")
include("markers/markers.jl")
include("variance.jl")
include("MCMCsamples.jl")
include("interface.jl")
include("Pi.jl")

export build_model,set_covariate,set_random,add_markers,get_pedigree
export showMME
export outputMCMCsamples,solve,runMCMC

end
