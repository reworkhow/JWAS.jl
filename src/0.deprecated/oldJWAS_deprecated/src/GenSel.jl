module GenSel

using(Distributions)

type GenSelOptions
    run          #BayeB,BayesC,BayesC0,BayesR
    seed         # set the seed for the random number generator
    chainLength  # number of iterations
    probFixed    # parameter "pi" the probability SNP effect is zero
    estimatePi   # "yes" or "no"
    estimateScale   # "yes" or "no"
    dfEffectVar  # hyper parameter (degrees of freedom) for locus effect variance
    nuRes        # hyper parameter (degrees of freedom) for residual variance
    varGenotypic # used to derive hyper parameter (scale) for locus effect variance
    varResidual  # used to derive hyper parameter (scale) for locus effect variance

    # belows are for BayesN
#    fittedSNPperWindow
#    windowWidth  # in mega-basepaire unit
#    markerMap    # marker_id, chrom_id, position

    #function GenSelOptions(run,seed,chainLength,probFixed,estimatePi,dfEffectVar,nuRes,varGenotypic,varResidual)#not necessary
    #    scaleVar = varGenotypic*(dfEffectVar-2)/dfEffectVar
    #    scaleRes = varResidual*(nuRes-2)/nuRes
    #    new(run,seed,chainLength,probFixed,estimatePi,dfEffectVar,nuRes,varGenotypic,varResidual,scaleVar,scaleRes)
    #end
end


export GenSelOptions #GenSel.options is not needed

include("Tools.jl")
include("Samplers.jl")
include("BayesC0.jl")
include("BayesCPi.jl")
include("BayesCPiDom.jl")
include("BayesB.jl")
include("BayesN.jl")
#include("BayesR.jl")



function runGenSel(myOptions::GenSelOptions,X,y,C,Rinv)
  srand(myOptions.seed)
  if myOptions.run=="BayesC0"
    BayesC0!(myOptions,X,y)
  elseif myOptions.run=="BayesB"
    BayesB!(myOptions,X,y,C,Rinv)
  elseif myOptions.run=="BayesCPi"
    BayesCPi!(myOptions,X,y,C,Rinv)
  elseif myOptions.run=="BayesR"
    #BayesR!(myOptions,X,y)
    BayesC0!(myOptions,X,y)
  elseif myOptions.run=="BayesCPiDom"
    BayesCPiDom!(myOptions,X,y,C,Rinv)
  elseif myOptions.run=="BayesN"
    BayesN!(myOptions,X,y,C,Rinv)
  end
end

export runGenSel
export get_dom_cov

end
