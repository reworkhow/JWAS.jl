# JWAS.jl

[![Build Status](https://travis-ci.org/reworkhow/JWAS.jl.svg?branch=master)](https://travis-ci.org/reworkhow/JWAS.jl)

JWAS.jl is an open-source software tool written in Julia for Bayesian multiple regression methods applied to genome-wide association studies and genomic prediction.

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.add("JWAS")`
* **Documentaion**: [available here](http://nbviewer.jupyter.org/github/reworkhow/JWAS.jl/tree/master/docs/index.ipynb)



### Structure of JWAS Module



```
JWAS.jl

├────── build_model
├────── set_covariate
├────── set_random
├────── get_pedigree
├────── add_markers
├────── outputMCMCsamples
├────── showMME
├────── solve
└────── runMCMC

├──────── PedModule.jl

├──────── Datasets.jl

├──────── SSBR.jl

├──────── misc.jl
           ├── get_breeding_values
           ├── get_additive_genetic_variances
           ├── get_heritability
           ├── get_correlations
           └── report


```

### Help

1. Show this README file in REPL or notebook using `?JWAS`
2. For help on a specific function above, type ? followed by its name, e.g. `?runMCMC` and press enter.
3. Run `Pkg.checkout("JWAS")` to get the newest unofficial JWAS. Run `Pkg.free("JWAS")` to go back to the offical one.

# SSBR

SSBR is a tool for single step Bayesian regression analyses.


####Quick-start

```Julia
using JWAS: Datasets,SSBR,misc

#data files from QTLDatasets package
pedfile    = Datasets.dataset("testSSBR","ped.txt")
genofile   = Datasets.dataset("testSSBR","genotype.txt")
phenofile  = Datasets.dataset("testSSBR","phenotype.txt")
fixedfile  = Datasets.dataset("testSSBR","fixed.txt")
Validation = Datasets.dataset("testSSBR","validation.txt")

#set up input parameters
input=InputParameters()
input.method       = "BayesC"
input.varGenotypic = 4.48
input.varResidual  = 6.72
input.probFixed    = 0.99
input.outFreq      = 10000


MCMCinfo(input)
#MCMC Information:
#seed                        314
#chainLength               50000
#method                   BayesC
#outFreq                    1000
#probFixed                 0.990
#varGenotypic              4.480
#varResidual               6.720
#estimateVariance           true
#estimatePi                false
#estimateScale             false
#dfEffectVar               4.000
#nuRes                     4.000
#nuGen                     4.000
#centering                 false


#run it
out=runSSBR(input,pedigree=pedfile,genotype=genofile,phenotype=phenofile,fixedfile=fixedfile);

#check accuracy
using DataFrames
df = readtable(Validation, eltypes =[UTF8String, Float64], separator = ' ',header=false,names=[:ID,:EBV]);
comp=join(out,df,on=:ID);
cor(comp[:EBV],comp[:EBV_1])

```
