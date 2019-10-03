# SSBR

SSBR is a tool for single step Bayesian regression analyses.


#### Quick-start

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
out=runSSBR(input,pedigree=pedfile,genotype=genofile,phenotype=phenofile,fixedfile=fixedfile); #return matrices, marker effects and ebv

#check accuracy
using DataFrames
df = readtable(Validation, eltypes =[String, Float64], separator = ' ',header=false,names=[:ID,:EBV]);
comp=join(out.ebv,df,on=:ID);
cor(comp[:EBV],comp[:EBV_1])
```

#### More

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Documentation**: [available here](https://github.com/QTL-rocks/SSBR.jl/wiki)
* **Authors**: [Hao Cheng](http://QTL.rocks), [Rohan Fernando](http://www.ans.iastate.edu/faculty/index.php?id=rohan)
