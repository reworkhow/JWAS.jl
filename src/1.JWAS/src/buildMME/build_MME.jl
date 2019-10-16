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
    build_model(model_equations::AbstractString,R=false; df::Float64=4.0)

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
function build_model(model_equations::AbstractString, R = 0.0; df = 4)
  if !isposdef(map(Float64,R)) && R != 0.0
    error("The covariance matrix is not positive definite.")
  end

  if !(typeof(model_equations)<:AbstractString) || model_equations==""
      error("Model equations are wrong.\n
      To find an example, type ?build_model and press enter.\n")
  end

  #e.g., ""y2 = A+B+A*B""
  modelVec   = [strip(i) for i in split(model_equations,[';','\n'],keepempty=false)]
  nModels    = size(modelVec,1)
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
  if R == 0.0
    R = (nModels==1 ? 0.0 : zeros(nModels,nModels))
  end
  return MME(nModels,modelVec,modelTerms,dict,lhsVec,map(Float64,R),Float64(df))
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
  trm.str = Array{AbstractString}(undef,nObs)
  trm.val = Array{Float64}(undef,nObs)

  for i = 1:trm.nFactors
    if trm.factors[i] != :intercept && any(ismissing,df[!,trm.factors[i]])
      error("Missing values are found in independent variables: ",trm.factors[i])
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
  trm.val = val
end

getFactor(str) = [strip(i) for i in split(str,"*")]

################################################################################
# make incidence matrix for each ModelTerm
#
################################################################################
function getX(trm::ModelTerm,mme::MME)
    #Row Index
    nObs  = length(trm.str)
    xi    = (trm.iModel-1)*nObs .+ collect(1:nObs)
    #Value
    xv    = trm.val
    #Column Index
    if trm.trmStr in mme.pedTrmVec
       #########################################################################
       #random polygenic effects,e.g."Animal","Animal*age"
       #column index needs to compromise numerator relationship matrix
       #########################################################################
       trm.names   = PedModule.getIDs(mme.ped)
       trm.nLevels = length(mme.ped.idMap)
       whichobs    = 1

       xj          = []
       for i in trm.str
         for animalstr in getFactor(i) #two ways:animal*age;age*animal
           if haskey(mme.ped.idMap, animalstr)   #"animal" ID not "age"
             xj = [xj; mme.ped.idMap[animalstr].seqID]
           elseif animalstr=="0" #founders "0" are not fitted effect (e.g, fitting maternal effects)
                                 #all non-founder animals in pedigree are effects
             xj = [xj; 1]        #put 1<=any interger<=nAnimal is okay)
             xv[whichobs]=0      #thus add one row of zeros in X
           end
         end
         whichobs    += 1
       end
       xj        = round.(Int64,xj)

       #some animal IDs in pedigree may be missing in data (df),ensure #columns = #animals in
       #pedigree by adding a column of zeros
       pedSize = length(mme.ped.idMap)
       xi      = [xi;1]          # adding a zero to
       xj      = [xj;pedSize]    # the last column in row 1
       xv      = [xv;0.0]
    else
       #########################################################################
       #other fixed or random effects
       #########################################################################
       #dict,trm.names  = mkDict(trm.str) #key: levels of variable; value: column index
       #trm.nLevels     = length(dict)
       #xj              = round.(Int64,[dict[i] for i in trm.str]) #column index
       res=Array{AbstractString,1}()
       for i in mme.rndTrmVec
           res=[res;i.term_array]
       end #make an array of all random effects (string)
       if trm.trmStr in res && mme.modelTermDict[trm.trmStr].names!=[]
           #imputational residual
           dict,thisnames  = mkDict(mme.modelTermDict[trm.trmStr].names)
           if thisnames!=mme.modelTermDict[trm.trmStr].names
               error("errors in SSBR")
           end
           if !issubset(filter(x->x≠"0",trm.str),thisnames)
             error("For trait ",trm.iTrait," some levels for ",trm.trmStr," in the phenotypic file are not found in levels for random effects ",
             trm.trmStr,". ","This may happen if the type is wrong, e.g, use of float instead of string.")
           end
       else
           dict,trm.names  = mkDict(trm.str)
           #delete "0"? for fixed effects missing
       end

       trm.nLevels     = length(dict)
       dict["0"]       = 1 #for missing data
       xj              = round.(Int64,[dict[i] for i in trm.str]) #column index
       xv[trm.str.=="0"] .= 0 #for missing data

       xi      = [xi;1]              # adding a zero to
       xj      = [xj;trm.nLevels]    # the last column in row 1
       xv      = [xv;0.0]
    end

    #ensure X has nObs*nModels rows
    nModels = size(mme.lhsVec,1)
    xi = [xi;nObs*nModels] #if (1,1) and (nObs*nModels,1) already exist
    xj = [xj;1]            # add 0 to X
    xv = [xv;0]

    #create X
    trm.X = sparse(xi,xj,xv)
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
    obsID = map(string,df[!,1])
    y   = recode(df[!,mme.lhsVec[1]], missing => 0.0)
    for i=2:size(mme.lhsVec,1)
      y   = [y; recode(df[!,mme.lhsVec[i]],missing=>0.0)]
    end
    ii    = 1:length(y)
    jj    = ones(length(y))
    vv    = y
    ySparse = sparse(ii,jj,vv)

    #make default covariance matrices if not provided
    set_default_priors_for_variance_components(mme,df)
    #Make lhs and rhs for MME
    mme.X       = X
    mme.ySparse = ySparse
    mme.obsID   = obsID

    if mme.nModels==1     #single-trait (lambda version)
      mme.mmeLhs = X'X
      mme.mmeRhs = X'ySparse
    elseif mme.nModels>1  #multi-trait
      Ri         = mkRi(mme,df) #handle missing phenotypes with ResVar
                                #make MME without variance estimation (constant)
                                #and residual imputation
      mme.mmeLhs = X'Ri*X
      mme.mmeRhs = X'Ri*ySparse
    end

    #Random effects parts in MME
    #Pedigree
    if mme.pedTrmVec != 0
        if mme.Ai == 0 #if no SSBR
            mme.Ai=PedModule.AInverse(mme.ped)
        end
      mme.GiOld = zeros(size(mme.GiOld))
      addA(mme::MME)
      mme.GiOld = copy(mme.GiNew)
    end

    #trick to enable addLambdas() 1st time
    #random_term.GiNew*mme.RNew - random_term.GiOld*mme.ROld
    for random_term in mme.rndTrmVec
      random_term.GiOld = zeros(size(random_term.GiOld))
    end
    addLambdas(mme)
    for random_term in mme.rndTrmVec
      random_term.GiOld = copy(random_term.GiNew)
    end

    dropzeros!(mme.mmeLhs)
    dropzeros!(mme.mmeRhs)
end

#set default covariance matrices for variance components if not provided
function set_default_priors_for_variance_components(mme,df)
  myvar     = [var(skipmissing((df[!,mme.lhsVec[i]]))) for i=1:size(mme.lhsVec,1)]
  phenovar  = diagm(0=>myvar)
  var_piece = phenovar/2

  #genetic variance or marker effect variance
  if mme.M!=0 && mme.M.G == false && mme.M.genetic_variance == false
    printstyled("Prior information for genomic variance is not provided and is generated from the data.\n",bold=false,color=:green)
    if mme.nModels==1
      mme.M.genetic_variance = var_piece[1,1]
    elseif mme.nModels>1
      mme.M.genetic_variance = var_piece
    end
  end
  #residual effects
  if mme.nModels==1 && isposdef(mme.RNew) == false #single-trait
    printstyled("Prior information for residual variance is not provided and is generated from the data.\n",bold=false,color=:green)
    mme.RNew = mme.ROld = var_piece[1,1]
  elseif mme.nModels>1 && isposdef(mme.R) == false #multi-trait
    printstyled("Prior information for residual variance is not provided and is generated from the data.\n",bold=false,color=:green)
    mme.R = var_piece
  end
  #polyginic effects
  if mme.pedTrmVec != 0 && isposdef(mme.Gi) == false
    printstyled("Prior information for polygenic effect variance is not provided and is generated from the data.\n",bold=false,color=:green)
    myvarout  = [split(i,":")[1] for i in mme.pedTrmVec]
    myvarin   = string.(mme.lhsVec)
    Zdesign   = mkmat_incidence_factor(myvarout,myvarin)
    G         = diagm(Zdesign*diag(var_piece))
    mme.Gi    = mme.GiOld = mme.GiNew = Symmetric(inv(G))
  end
  #other random effects
  if length(mme.rndTrmVec) != 0
    for randomEffect in mme.rndTrmVec
      if isposdef(randomEffect.Gi) == false
        printstyled("Prior information for random effect variance is not provided and is generated from the data.\n",bold=false,color=:green)
        myvarout  = [split(i,":")[1] for i in randomEffect.term_array]
        myvarin   = string.(mme.lhsVec)
        Zdesign   = mkmat_incidence_factor(myvarout,myvarin)
        G         = diagm(Zdesign*diag(var_piece))

        randomEffect.Gi = randomEffect.GiOld = randomEffect.GiNew = Symmetric(inv(G))
        randomEffect.scale = G*(randomEffect.df-length(randomEffect.term_array)-1)
      end
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
