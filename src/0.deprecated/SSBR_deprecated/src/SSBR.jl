module SSBR

using ..misc
using ..PedModule
using DataFrames
using SparseArrays
using Printf

include("SSBR_types.jl")
include("getMatrices.jl")
include("imputationResidual/sampleEpsilon.jl")
include("ssBayesianAlphabet/ssBayesC0.jl")
include("ssBayesianAlphabet/ssBayesC0_constantvariance.jl")#for test
include("ssBayesianAlphabet/ssBayesC_constantvariance.jl")#for test
include("ssBayesianAlphabet/ssBayesC.jl")
include("ssBayesianAlphabet/ssBayesB.jl")
include("run.jl")

#include("ssMME.jl")

export runSSBR
#export ssBayesC0_constantvariance
#export ssBayesB,ssBayesC,ssBayesC0
#export ssMME
#export PBLUP


end
