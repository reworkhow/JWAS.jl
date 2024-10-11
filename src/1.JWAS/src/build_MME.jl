"""
    mkDict(a::Vector{T}) where T <: Any

    Get column index in the incidence matrix for each level of a factor (categorical variable) 
    input:  a=["a1","a4","a1","a2"] 
    output: d=Dict("a2" => 3, "a1" => 1, "a4" => 2), level_names=["a1","a4","a2"]
    
    note: enumerate(level_names) gives a list of tuples (index, element), reverse() to reverse (index,element) to (element,index)
"""
function mkDict(a::Vector{T}) where T <: Any
    # Create a vector of unique levels from a
    level_names = String.(unique(a)) #e.g., ["a1","a2","a1"] -> ["a1","a2"]; apply String() to avoid the String3 type in enumerate()
    # Create a dictionary mapping each element to its index in level_names
    d = Dict(map(reverse, enumerate(level_names))) #e.g., Dict("a1"=>1,"a2"=>2), level_names=["a1","a2"]
    return d, level_names
end

"""
    build_model(model_equations::AbstractString,R=false; df::AbstractFloat=4.0, estimate_variance=true)

* Build a model from **model equations** with the residual variance **R**. In Bayesian analysis, **R**
  is the mean for the prior assigned for the residual variance with degree of freedom **df**, defaulting
  to 4.0. If **R** is not provided, a value is calculated from responses (phenotypes).
* By default, all variabels in model_equations are factors (categorical) and fixed. Set variables
  to be covariates (continuous) or random using functions `set_covariate()` or `set_random()`.
* The argument `estimate_variance` indicates whether to estimate the residual variance; `estimate_variance=true` is the default.

```julia
#single-trait
model_equations = "BW = intercept + age + sex"
R               = 6.72
models          = build_model(model_equations,R);

#multi-trait
model_equations = "BW = intercept + age + sex
                   CW = intercept + litter";
R               = [6.72   24.84
                   24.84  708.41]
models          = build_model(model_equations,R);
```
"""
function build_model(model_equations::AbstractString, 
                     ## residual variance: 
                     R = false; df = 4.0, 
                     estimate_variance=true, estimate_scale=false, 
                     constraint=false, #for multi-trait only, constraint=true means no residual covariance among traits
                     ## censored, categorical traits:
                     censored_trait = false, categorical_trait = false)

    if R != false && !isposdef(map(AbstractFloat,R))
      error("The covariance matrix is not positive definite.")
    end
    if !(typeof(model_equations)<:AbstractString) || model_equations==""
      error("Model equations are wrong.\n
      To find an example, type ?build_model and press enter.\n")
    end
    if estimate_scale != false
      error("estimate scale for residual variance is not supported now.")
    end

    ############################################################################
    # All model terms (will be added to MME)
    ############################################################################
    #e.g., ""y2 = A+B+A*B""
    modelVec   = [strip(i) for i in split(model_equations,[';','\n'],keepempty=false)]
    nModels    = size(modelVec,1)
    if R != false && size(R,1) != nModels
      error("The residual covariance matrix is not a ",nModels," by ",nModels," matrix.")
    end
    lhsVec     = Symbol[]    #:y, phenotypes
    modelTerms = ModelTerm[] #initialization of an array of ModelTerm outside for loop
    dict       = Dict{AbstractString,ModelTerm}()
    for (m,model) = enumerate(modelVec)
      lhsRhs = split(model,"=")                  #"y2","A+B+A*B"
      lhs    = strip(lhsRhs[1])                  #"y2"
      lhsVec = [lhsVec;Symbol(lhs)]              #:y2
      rhsVec = split(strip(lhsRhs[2]),"+")       #"A","B","A*B"
      mTrms  = [ModelTerm(strip(trmStr),m,lhs) for trmStr in rhsVec]
      modelTerms  = [modelTerms;mTrms]           #a vector of ModelTerm
    end
    for trm in modelTerms          #make a dict for model terms
      dict[trm.trmStr] = trm
    end

    ############################################################################
    # Genotypes (will be added to MME)
    ############################################################################
    genotypes = []
    whichterm = 1
    for term in modelTerms
      term_symbol = Symbol(split(term.trmStr,":")[end])
      if isdefined(Main,term_symbol) #@isdefined can be used to test whether a local variable or object field is defined
        if typeof(getfield(Main,term_symbol)) == Genotypes
          term.random_type = "genotypes"
          genotypei = getfield(Main,term_symbol)
          genotypei.name = string(term_symbol)
          trait_names=[term.iTrait]
          if genotypei.name ∉ map(x->x.name, genotypes) #only save unique genotype
            genotypei.ntraits = nModels
            genotypei.trait_names = string.(lhsVec)
            if nModels != 1
              genotypei.G.df = genotypei.G.df + nModels
            end
            if (genotypei.G.val != false || genotypei.genetic_variance.val != false)
              if size(genotypei.G.val,1) != nModels && size(genotypei.genetic_variance.val,1) != nModels
                error("The genomic covariance matrix is not a ",nModels," by ",nModels," matrix.")
              end
            end
            push!(genotypes,genotypei)
          end
        end
      end
    end

  #create mme with genotypes
  filter!(x->x.random_type != "genotypes",modelTerms) #remove "genotypes" from modelTerms
  filter!(x->x[2].random_type != "genotypes",dict)    #remove "genotypes" from dict

  #set scale and df for residual variance
  if nModels == 1
    scale_R = R*(df - 2)/df
    df_R    = df 
  else
    scale_R = R*(df - 1)
    df_R    = df + nModels
  end
  
  #initialize mme
  mme = MME(nModels,modelVec,modelTerms,dict,lhsVec, 
            Variance(R==false ? R : Float32.(R), #val
                     Float32(df_R),              #df
                     R==false ? R : scale_R,     #scale
                     estimate_variance, estimate_scale, constraint))
  if length(genotypes) != 0
    mme.M = genotypes #add genotypes into mme
  end

  #setup traits_type (by default is "continuous")
  for t in 1:mme.nModels
    if string(mme.lhsVec[t]) ∈ censored_trait
      mme.traits_type[t]="censored"
    elseif string(mme.lhsVec[t]) ∈ categorical_trait
      mme.traits_type[t]="categorical"
    end
  end

  return mme
end

"""
    set_covariate(model::MME,variables::AbstractString...)

* set **variables** as covariates; **model** is the output of function `build_model()`.

```julia
#After running build_model, variabels age and year can be set to be covariates as
set_covariate(model,"age","year")
#or
set_covariate(model,"age year")
```
"""
function set_covariate(mme::MME,covStr::AbstractString...)
  covVec=[]
  for i in covStr
    covVec = [covVec;strip.(split(i," ",keepempty=false))]
  end
  mme.covVec = [mme.covVec;[Symbol(i) for i in covVec]]
end

################################################################################
#Get all data from data files (in DataFrame) based on each ModelTerm
#Fill up str and val for each ModelTerm
################################################################################

function getData(trm::ModelTerm,df::DataFrame,mme::MME) #ModelTerm("1:A*B")
  nObs    = size(df,1)
  if trm.factors[1] == :intercept #for intercept
    str = fill("intercept",nObs)
    val = fill(1.0,nObs)
  else                            #for ModelTerm e.g. "1:A*B" (or "1:A")
    myDf = df[!,trm.factors]                        #:A,:B
    if trm.factors[1] in mme.covVec                 #if A is a covariate
      if !(typeof(df[!,trm.factors[1]][1]) <: Number)
        error("$(trm.factors[1]) is fitted as a covariate (continuous variables). The data type should be numbers.")
      end
      str = fill(string(trm.factors[1]),nObs)       #["A","A",...]
      val = df[!,trm.factors[1]]                    #df[:A]
    else                                              #if A is a factor (animal or maternal effects)
      str = [string(i) for i in df[!,trm.factors[1]]] #["A1","A3","A2","A3",...]
      val = fill(1.0,nObs)                            #1.0,1.0...
    end

    #for ModelTerm object e.g. "A*B" whose nFactors>1
    for i=2:trm.nFactors
      if trm.factors[i] in mme.covVec
        #["A * B","A * B",...] or ["A1 * B","A2 * B",...]
        if !(typeof(df[!,trm.factors[i]][1]) <: Number)
          error("$(trm.factors[i]) is fitted as a covariate (continuous variables). The data type should be numbers.")
        end
        str = str .* fill(" * "*string(trm.factors[i]),nObs)
        val = val .* df[!,trm.factors[i]]
      else
        #["A * B1","A * B2",...] or ["A1 * B1","A2 * B2",...]
        str = str .* fill(" * ",nObs) .* [string(j) for j in df[!,trm.factors[i]]]
        val = val .* fill(1.0,nObs)
      end
    end
  end
  trm.data = str
  val=convert(Array,val)
  val = coalesce.(val, 0.0) #replace missing with 0.0
  trm.val = ((mme.MCMCinfo == false || mme.MCMCinfo.double_precision) ? Float64.(val) : Float32.(val))
end

getFactor(str) = [strip(i) for i in split(str,"*")]

################################################################################
# make incidence matrix for each ModelTerm
#
################################################################################
function getX(trm::ModelTerm,mme::MME)
    #1. Row Index
    nObs  = length(trm.data)
    xi    = (trm.iModel-1)*nObs .+ collect(1:nObs)
    #2. Value
    xv    = trm.val
    #3. Column Index
    #replace missing values (e.g., "missing","missing*a1") with "missing"
    for i in 1:length(trm.data)  #data
      if "missing" in getFactor(trm.data[i]) #replace "A*missing" with "missing"
          trm.data[i] = "missing"
      end
    end
    #3.1. fixed effects or iid random effects (names: "a1","a1*b1")
    if trm.random_type == "I" || trm.random_type == "fixed"
      dict,trm.names = mkDict(filter(x->x≠"missing",trm.data))
      trm.nLevels    = length(dict)
      data           = trm.data
    end
    if trm.random_type == "V" || trm.random_type == "A"
      #########################################################################
      #random polygenic effects,e.g."Animal","Animal*age"
      #column index needs to compromise numerator relationship matrix
      #########################################################################
       dict,trm.names  = mkDict(mme.modelTermDict[trm.trmStr].names) #key: levels of variable; value: column index
       trm.nLevels  = length(dict) #before add key "missing"
       #3.2. Levels for each observation
       #Get the random effect in interactions
       data=[]
       for i in trm.data
         for factorstr in getFactor(i) #two ways:animal*age;age*animal
           if factorstr in trm.names || factorstr == "missing"  #"animal" ID not "age"
             data = [data;factorstr]
           end
         end
       end
       if length(data) < length(trm.data)
         error("For trait ",trm.iTrait," some levels for ",trm.trmStr," in the phenotypic file are not found in levels for random effects ",
         trm.trmStr,". ","This may happen if missing values are not considered in missingstrings.")
       elseif length(data) > length(trm.data)
         error("Same level names are found for the two terms in the interaction.")
       end
    end
    #4. missing values in random effects
    #e.g, founders in pedigree are "missing" ("0") for materal effects
    #e.g, ϵ in single-step SNP-BLUP are missings for genotyped individuals
    #Thus, add one row of zeros in the design matrix for corresponding observation
    dict["missing"]          = 1 #column index for missing data (any positive integer < nLevels)
    xv[data.=="missing"]     .= 0 #values       for missing data

    #5.column index
    xj                 = round.(Int64,[dict[i] for i in data]) #column index

    #create the design matrix
    #ensure size of X is nObs*nModels X nLevels
    nModels = size(mme.lhsVec,1)
    #create X
    trm.X = sparse(xi,xj,xv,nObs*nModels,trm.nLevels)
    dropzeros!(trm.X)
    trm.startPos = mme.mmePos
    mme.mmePos  += trm.nLevels
end

"""
Construct mixed model equations with

incidence matrix: X      ;
response        : ySparse;
left-hand side  : mmeLhs ;
right-hand side : mmeLhs ;
"""
function getMME(mme::MME, df::DataFrame)
    if mme.mmeLhs != false
      error("Please build your model again using the function build_model().")
    end
    #Heterogeneous residuals
    if mme.MCMCinfo != false && mme.MCMCinfo.heterogeneous_residuals == true
        invweights = 1 ./ convert(Array,df[!,Symbol("weights")])
    else
        invweights = ones(size(df,1))
    end
    mme.invweights = (mme.MCMCinfo == false || mme.MCMCinfo.double_precision ? Float64.(invweights) : Float32.(invweights))
    #Make incidence matrices X for each term
    for trm in mme.modelTerms
      if trm.X == false
        getData(trm,df,mme)
        getX(trm,mme)
      end
    end
    #concatenate all terms
    X   = mme.modelTerms[1].X
    for i=2:length(mme.modelTerms)
       X = [X mme.modelTerms[i].X]
    end


    y   = coalesce.(df[!,mme.lhsVec[1]], 0.0) #replace missing values with 0
    for i=2:size(mme.lhsVec,1)
      y   = [y; coalesce.(df[!,mme.lhsVec[i]], 0.0)]
    end
    ii      = 1:length(y)
    jj      = ones(length(y))
    vv      = ((mme.MCMCinfo == false || mme.MCMCinfo.double_precision) ? Float64.(y) : Float32.(y))
    ySparse = sparse(ii,jj,vv)

    #Make lhs and rhs for MME
    mme.X       = X
    mme.ySparse = ySparse

    if mme.nModels==1     #single-trait (lambda version)
        mme.mmeLhs = X'*Diagonal(mme.invweights)*X
        mme.mmeRhs = X'*Diagonal(mme.invweights)*ySparse
    elseif mme.nModels>1  #multi-trait
    #Ri of size (nind x ntraits)^2, the inverse of the variance of the random
    #vector of errors, is generated based on the phenotype missing patterns,
    #such that no imputation of missing phenotypes is required.
    #mixed model equations is obtained below for multi-trait PBLUP
    #with known residual covariance matrix and missing phenotypes.
      # if mme.MCMCinfo.mega_trait == true  #multiple single trait
      if mme.R.constraint == true #tj: Hao, please confirm! now we do not have mege_trait option. We only have constraint option for variances
        Ri = Diagonal(repeat(mme.invweights,mme.nModels))
      else  #multi-trait
        Ri = mkRi(mme,df,mme.invweights)
      end
      mme.mmeLhs = X'Ri*X
      mme.mmeRhs = X'Ri*ySparse
    end

    #Random effects parts in MME
    if mme.nModels == 1
      #random_term.GiNew*mme.R.val - random_term.GiOld*mme.ROld
      for random_term in mme.rndTrmVec #trick
        random_term.GiOld.val = zero(random_term.GiOld.val)
      end
      addVinv(mme)
      for random_term in mme.rndTrmVec #trick
        random_term.GiOld.val = copy(random_term.GiNew.val)
      end
    else
      addVinv(mme)
    end

    dropzeros!(mme.mmeLhs)
    dropzeros!(mme.mmeRhs)

    #No phenotypic data for some levels of a factor in multi-trait analysis
    #e.g., y3:x3:f in https://github.com/reworkhow/JWAS.jl/blob/
    #a6f4595796b70811c0b745af525e7e0a822bb954/src/5.Datasets/data/example/phenotypes.txt
    for i in size(mme.mmeLhs,1)
      if mme.mmeLhs[i,i] == 0.0
        error("No phenotypic data for ",getNames(mme)[i])
      end
    end
end

#Get names for variables in Mixed Model Equations in order
function getNames(mme::MME)
    names = Array{AbstractString}(undef,0)
    for trm in mme.modelTerms
        for name in trm.names
            push!(names,trm.trmStr*":"*name)
        end
    end
    return names
end

function getNames(trm::ModelTerm)
    names = Array{AbstractString}(undef,0)
    for name in trm.names
        push!(names,trm.trmStr*":"*name)
    end
    return names
end
