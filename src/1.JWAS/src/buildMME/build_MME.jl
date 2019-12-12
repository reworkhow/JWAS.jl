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
  is the mean for the prior assigned for the residual variance with degree of freedom **df** defaulting
  to 4.0. If **R** is not provided, a value is calculate from responses (phenotypes).
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
function build_model(model_equations::AbstractString, R = false; df = 4.0)
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
  return MME(nModels,modelVec,modelTerms,dict,lhsVec,R == false ? R : Float32.(R),Float32(df))
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
  for i = 1:trm.nFactors
    if trm.factors[i] != :intercept && any(ismissing,df[!,trm.factors[i]])
      printstyled("Missing values are found in independent variables: ",trm.factors[i],".\n",bold=false,color=:red)
    end
  end

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
  trm.str = str
  trm.val = (mme.MCMCinfo.double_precision ? Float64.(val) : Float32.(val))
end

getFactor(str) = [strip(i) for i in split(str,"*")]

################################################################################
# make incidence matrix for each ModelTerm
#
################################################################################
function getX(trm::ModelTerm,mme::MME)
    #1. Row Index
    nObs  = length(trm.str)
    xi    = (trm.iModel-1)*nObs .+ collect(1:nObs)
    #2. Value
    xv    = trm.val
    #3. Column Index
    #3.1. position fo all levels
    if trm.random_type == "V" || trm.random_type == "A"
      #########################################################################
      #random polygenic effects,e.g."Animal","Animal*age"
      #column index needs to compromise numerator relationship matrix
      #########################################################################
       dict,trm.names  = mkDict(mme.modelTermDict[trm.trmStr].names) #key: levels of variable; value: column index
       if trm.names!=mme.modelTermDict[trm.trmStr].names
           error("The order of names is changed!")
       end
       if trm.nFactors == 1 && !issubset(filter(x->x≠"missing",trm.str),trm.names)
         error("For trait ",trm.iTrait," some levels for ",trm.trmStr," in the phenotypic file are not found in levels for random effects ",
         trm.trmStr,". ","This may happen if the type is wrong, e.g, use of float instead of string.")
       end
    else #fixed or iid random effects (also works with interactions)
      dict,trm.names  = mkDict(trm.str)
    end
    trm.nLevels  = length(dict) #before add key "missing"
    #3.2. Levels for each observation
    #Get the random effect in interactions
    if (trm.random_type == "V" || trm.random_type == "A") && trm.nFactors != 1
      str=[]
      for i in trm.str  #data
        for factorstr in getFactor(i) #two ways:animal*age;age*animal
          if factorstr in trm.names || factorstr == "missing"  #"animal" ID not "age"
            str = [str;factorstr]
          end
        end
      end
      if length(str) != length(trm.str)
        error("Same level names are found in the two factors in the interaction.")
      end
      trm.str=str
    end
    #4. missing values in random effects
    #e.g, founders in pedigree are "missing" ("0") for materal effects
    #e.g, ϵ in single-step SNP-BLUP are missings for genotyped individuals
    #Thus, add one row of zeros in the design matrix for corresponding observation
    dict["missing"]          = 1 #column index for missing data (any positive integer < nLevels)
    xv[trm.str.=="missing"] .= 0 #values       for missing data

    #5.column index
    xj                 = round.(Int64,[dict[i] for i in trm.str]) #column index

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
    if mme.mmePos != 1
      error("Please build your model again using the function build_model().")
    end

    #Make incidence matrices X for each term
    for trm in mme.modelTerms
      getData(trm,df,mme)
      getX(trm,mme)
    end
    #concatenate all terms
    X   = mme.modelTerms[1].X
    for i=2:length(mme.modelTerms)
       X = [X mme.modelTerms[i].X]
    end

    #Make response vector (y)
    y   = recode(df[!,mme.lhsVec[1]], missing => 0.0)
    for i=2:size(mme.lhsVec,1)
      y   = [y; recode(df[!,mme.lhsVec[i]],missing=>0.0)]
    end
    ii      = 1:length(y)
    jj      = ones(length(y))
    vv      = (mme.MCMCinfo.double_precision ? Float64.(y) : Float32.(y))
    ySparse = sparse(ii,jj,vv)

    #Make lhs and rhs for MME
    mme.X       = X
    mme.ySparse = ySparse

    if mme.nModels==1     #single-trait (lambda version)
      if mme.invweights == false
        mme.mmeLhs = X'X
        mme.mmeRhs = X'ySparse
      else
        mme.mmeLhs = X'*Diagonal(mme.invweights)*X
        mme.mmeRhs = X'*Diagonal(mme.invweights)*ySparse
      end
    elseif mme.nModels>1  #multi-trait
      Ri         = mkRi(mme,df) #handle missing phenotypes with ResVar
                                #make MME without variance estimation (constant)
                                #and residual imputation
      mme.mmeLhs = X'Ri*X
      mme.mmeRhs = X'Ri*ySparse
    end

    #Random effects parts in MME
    #trick to enable addLambdas() 1st time ???
    #random_term.GiNew*mme.RNew - random_term.GiOld*mme.ROld
    for random_term in mme.rndTrmVec
      random_term.GiOld = zero(random_term.GiOld)
    end
    addVinv(mme)
    for random_term in mme.rndTrmVec
      random_term.GiOld = copy(random_term.GiNew)
    end

    dropzeros!(mme.mmeLhs)
    dropzeros!(mme.mmeRhs)
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
            push!(names,trm.trmStr*" : "*name)
        end
    end
    return names
end

function getNames(trm::ModelTerm)
    names = Array{AbstractString}(undef,0)
    for name in trm.names
        push!(names,trm.trmStr*" : "*name)
    end
    return names
end
