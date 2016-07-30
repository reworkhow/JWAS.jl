module QTL

using DataFrames
using Distributions
using PyPlot

#=
include("genotypes.jl")
include("fixed_effects.jl")
include("QTL_types.jl")
include("tools.jl")
include("files.jl")
include("QualityControl/QualityControl.jl")
include("Gibbs.jl")
include("BayesianAlphabet/BayesB.jl")
include("BayesianAlphabet/BayesC.jl")
=#
include("using_marker_samples.jl")
include("stat.jl")

end
