using Distributions,Printf,Random
using DelimitedFiles
using InteractiveUtils #for versioninfo
using DataFrames,CSV
using SparseArrays
using LinearAlgebra
using ProgressMeter
using .PedModule
using .misc

#Models
include("buildMME/types.jl")
include("buildMME/build_MME.jl")
include("buildMME/random_effects.jl")
include("buildMME/residual.jl")
include("buildMME/sample_variances.jl")
include("buildMME/solver.jl")
include("buildMME/model_info.jl")

#Markov chain Monte Carlo
include("MCMC/runMCMC.jl")
include("MCMC/outputMCMCsamples.jl")
include("MCMC/precheck.jl")
include("MCMC/MCMC_BayesC.jl")
include("MCMC/MCMC_GBLUP.jl")
include("MCMC/MT_MCMC_BayesC.jl")
include("MCMC/MT_PBLUP_constvare.jl")
include("MCMC/output.jl")

#Genomic Markers
include("markers/tools4genotypes.jl")
include("markers/readgenotypes.jl")
include("markers/BayesianAlphabet/BayesC0.jl")
include("markers/BayesianAlphabet/BayesC.jl")
include("markers/BayesianAlphabet/BayesB.jl")
include("markers/BayesianAlphabet/MTBayesC.jl")
include("markers/BayesianAlphabet/MTBayesCC.jl")
include("markers/BayesianAlphabet/MTBayesB.jl")
include("markers/BayesianAlphabet/MTBayesC0L.jl")
include("markers/Pi.jl")

#Incomplete Genomic Data (Single-step Methods)
include("SSBR/SSBR.jl")

#Structure Equation Models
include("StructureEquationModel/SEM.jl")

#MISC
include("misc/misc.jl")
include("pipeline/pipeline.jl")

export build_model,set_covariate,set_random,add_genotypes
export outputMCMCsamples,outputEBV,getEBV
export solve,runMCMC
export showMME,getinfo
#Pedmodule
export get_pedigree
#misc
export GWAS,QC
export get_additive_genetic_variances,get_breeding_values
export get_correlations,get_heritability
#others
export adjust_phenotypes,LOOCV
