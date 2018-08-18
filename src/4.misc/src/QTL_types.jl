using Printf 
mutable struct InputParameters
  seed::Int64            # seed
  method::AbstractString # BaysABC
  chainLength::Int64     # number of iterations
  outFreq::Int64         # frequency to print out MCMC status
  probFixed::Float64     # parameter "pi" the probability SNP effect is zero
  varGenotypic::Float64  # used to derive hyper parameter (scale) for locus effect variance
  varResidual::Float64   # used to derive hyper parameter (scale) for locus effect variance
  estimateVariance::Bool # true or false
  estimatePi::Bool       # true or false
  estimateScale::Bool    # true or false
  dfEffectVar::Float64   # hyper parameter (degrees of freedom) for locus effect variance
  nuRes::Float64         # hyper parameter (degrees of freedom) for residual variance
  nuGen::Float64         # hyper parameter (degree of freedom) for genetic variance (for ϵ)
  centering::Bool        # center real genotype or not
end

InputParameters()=InputParameters(314,"BayesC",50000,1000,0.95,1.0,1.0,true,false,false,4,4,4,false)

function MCMCinfo(input::InputParameters)
    println("MCMC Information:")
    @printf("%-20s %10s\n","seed",input.seed)
    @printf("%-20s %10s\n","chainLength",input.chainLength)
    @printf("%-20s %10s\n","method",input.method)
    @printf("%-20s %10d\n","outFreq",input.outFreq)
    @printf("%-20s %10.3f\n","probFixed",input.probFixed)
    @printf("%-20s %10.3f\n","varGenotypic",input.varGenotypic)
    @printf("%-20s %10.3f\n","varResidual",input.varResidual)
    @printf("%-20s %10s\n","estimateVariance",input.estimateVariance)
    @printf("%-20s %10s\n","estimatePi",input.estimatePi)
    @printf("%-20s %10s\n","estimateScale",input.estimateScale)
    @printf("%-20s %10.3f\n","dfEffectVar",input.dfEffectVar)
    @printf("%-20s %10.3f\n","nuRes",input.nuRes)
    @printf("%-20s %10.3f\n","nuGen",input.nuGen)
    @printf("%-20s %10s\n","centering",input.centering)
end

mutable struct GibbsMats #in src/1.JWAS/markers/tools
    X::Array{Float64,2}
    nrows::Int64
    ncols::Int64
    xArray::Array{Array{Float64,1},1}
    xpx::Array{Float64,1}
    function GibbsMats(X::Array{Float64,2}) ###More
        nrows,ncols = size(X)
        xArray = get_column_ref(X)
        XpX = getXpRinvX(X)
        new(X,nrows,ncols,xArray,XpX)
    end
end

mutable struct Current
    varGenotypic::Float64      #genotypic variance
    varResidual::Float64       #residual variance
    varEffect::Float64  # common marker variance
    scaleVar::Float64   # scale factor for locus effects
    scaleRes::Float64   # scale factor for residual varianc
    scaleGen::Float64
    fixed_effects       # sample of fixed effects β
    α                   # sample of partial marker effects unconditional on δ
    δ                   # inclusion indicator for marker effects
    u                   # sample of marker effects,u=α·δ
    π                   # probFixed
    locusEffectVar      #locus-specific variance

    iter
    nLoci               #count number of markers in the model
    yCorr

    imputation_residual #residual in SSBR

    function Current(input::InputParameters,geno::Genotypes,
                     fixed::FixedMatrix,y::Array{Float64,1})
        varGenotypic= input.varGenotypic
        varResidual = input.varResidual
        π           = input.probFixed
        dfEffectVar = input.dfEffectVar
        nuRes       = input.nuRes
        nuGen       = input.nuGen
        varGenotypic= input.varGenotypic
        sum2pq      = geno.sum2pq
        nMarkers    = geno.nMarkers
        nFixedEffects  = size(fixed.C,2)

        varEffect      = varGenotypic/((1-π)*sum2pq)
        locusEffectVar = fill(varEffect,nMarkers) #add if statement later
        scaleVar       = varEffect*(dfEffectVar-2)/dfEffectVar # scale factor for locus effects
        scaleRes       = varResidual*(nuRes-2)/nuRes                   # scale factor for residual varianc
        scaleGen       = varGenotypic*(nuGen-2)/nuGen                   # scale factor for residual varianc
        β          = zeros(nFixedEffects)  # sample of fixed effects
        α          = zeros(nMarkers)       # sample of partial marker effects unconditional on δ
        δ          = zeros(nMarkers)       # inclusion indicator for marker effects
        u          = zeros(nMarkers)       # sample of marker effects

        iter       = 0
        nLoci      = 0
        yCorr      = copy(y)
        ϵ          = zeros(1)

        new(varGenotypic,varResidual,varEffect,scaleVar,scaleRes,scaleGen,
            β,α,δ,u,π,locusEffectVar,iter,nLoci,yCorr,ϵ)
    end
end

mutable struct Output
    mean_fixed_effects::Array{Float64,1}
    meanMarkerEffects::Array{Float64,1}
    modelFreq::Array{Float64,1}
    resVar::Array{Float64,1}
    genVar::Array{Float64,1}
    pi::Array{Float64,1}
    scale::Array{Float64,1}

    mean_imputation_residual::Array{Float64,1}

    function Output(input::InputParameters,geno::Genotypes,fixed::FixedMatrix)
      chainLength   = input.chainLength
      nFixedEffects = size(fixed.C,2)#easy bug
      nMarkers      = geno.nMarkers
      new(zeros(nFixedEffects),
          zeros(nMarkers),
          zeros(nMarkers),
          zeros(chainLength),
          zeros(chainLength),
          zeros(chainLength),
          zeros(chainLength),
          zeros(1))
    end
end

export InputParameters
#export Output
#export GibbsMats
#export Current
export MCMCinfo
