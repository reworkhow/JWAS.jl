module misc

using DataFrames,CSV,DelimitedFiles
using Distributions
using ProgressMeter

include("genotypes.jl")
include("fixed_effects.jl")
include("QTL_types.jl")
include("tools.jl")
include("files.jl")
include("qc/qc.jl")
include("Gibbs.jl")
include("BayesianAlphabet/BayesB.jl")
include("BayesianAlphabet/BayesC.jl")
include("using_marker_samples.jl")
include("GWAS.jl")

end
