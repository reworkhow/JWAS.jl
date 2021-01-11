#return column index for each level of variables in incidence matrix X
#e.g. "A1"=>1,"A2"=>2
function mkDict(a)
  aUnique = unique(a)
  d = Dict()
  names = Array{Any}(undef,size(aUnique,1))
  for (i,s) in enumerate(aUnique)
    names[i] = s
    d[s] = i
  end
  return d,names
end

"""
    build_model(model_equations::AbstractString,R=false; df::AbstractFloat=4.0)

* Build a model from **model equations** with the residual variance **R**. In Bayesian analysis, **R**
  is the mean for the prior assigned for the residual variance with degree of freedom **df**, defaulting
  to 4.0. If **R** is not provided, a value is calculated from responses (phenotypes).
* By default, all variabels in model_equations are factors (categorical) and fixed. Set variables
  to be covariates (continuous) or random using functions `set_covariate()` or `set_random()`.

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
function build_model(model_equations::AbstractString, R = false; df = 4.0,
                     num_latent_traits = false, nonlinear_function = false) #nonlinear_function(x1,x2) = x1+x2
  if num_latent_traits != false
    lhs, rhs = strip.(split(model_equations,"="))
    model_equations = ""
    for i = 1:num_latent_traits
      model_equations = model_equations*lhs*string(i)*"="*rhs*";"
    end
    model_equations = model_equations[1:(end-1)]
  end

  if R != false && !isposdef(map(AbstractFloat,R))
    error("The covariance matrix is not positive definite.")
  end

  if !(typeof(model_equations)<:AbstractString) || model_equations==""
      error("Model equations are wrong.\n
      To find an example, type ?build_model and press enter.\n")
  end

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

  #add genotypes to mme
  genotypes = []
  whichterm = 1
  for term in modelTerms
    term_symbol = Symbol(split(term.trmStr,":")[end])
    traiti      = term.iModel
    if isdefined(Main,term_symbol) #@isdefined can be usde to tests whether a local variable or object field is defined
      if typeof(getfield(Main,term_symbol)) == Genotypes
        term.random_type = "genotypes"
        if traiti == 1 #same genos are required in all traits
          genotypei = getfield(Main,term_symbol)
          genotypei.name = string(term_symbol)
          genotypei.ntraits = nModels
          if nModels != 1
            genotypei.df = genotypei.df + nModels
          end
          if genotypei.G != false || genotypei.genetic_variance != false
            if size(genotypei.G,1) != nModels && size(genotypei.genetic_variance,1) != nModels
              error("The genomic covariance matrix is not a ",nModels," by ",nModels," matrix.")
            end
          end
          push!(genotypes,genotypei)
        end
      end
    end
  end
  #crear mme with genotypes
  filter!(x->x.random_type != "genotypes",modelTerms)
  mme = MME(nModels,modelVec,modelTerms,dict,lhsVec,R == false ? R : Float32.(R),Float32(df))
  if length(genotypes) != 0
    mme.M = genotypes
  end

  #laten traits
  if num_latent_traits != false
    mme.latent_traits = true
    if nonlinear_function != false
      mme.nonlinear_function = nonlinear_function
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
  DataFrames.recode!(val, missing => 0.0)
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

    #Make response vector (y)
    y   = DataFrames.recode(df[!,mme.lhsVec[1]], missing => 0.0)
    for i=2:size(mme.lhsVec,1)
      y   = [y; DataFrames.recode(df[!,mme.lhsVec[i]],missing=>0.0)]
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
      Ri         = mkRi(mme,df,mme.invweights)
      mme.mmeLhs = X'Ri*X
      mme.mmeRhs = X'Ri*ySparse
    end

    #Random effects parts in MME
    if mme.nModels == 1
      #random_term.GiNew*mme.R - random_term.GiOld*mme.ROld
      for random_term in mme.rndTrmVec #trick
        random_term.GiOld = zero(random_term.GiOld)
      end
      addVinv(mme)
      for random_term in mme.rndTrmVec #trick
        random_term.GiOld = copy(random_term.GiNew)
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

################################################################################
#Get left-hand side and right-hand side of the mixed model equation (no markers)
################################################################################
"""
    showMME(mme::MME,df::DataFrame)

* Show left-hand side and right-hand side of mixed model equations (no markers).
"""
function showMME(mme::MME,df::DataFrame)
   if size(mme.mmeRhs)==()
     getMME(mme,df)
   end
   return [getNames(mme) mme.mmeLhs],[getNames(mme) mme.mmeRhs]
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
