function mkDict(a)
    aUnique = unique(a)
    d = Dict()
    names = Array(Any,size(aUnique,1))
    for (i,s) in enumerate(aUnique)
        names[i] = s
        d[s] = i
    end
    return d,names
end

"""
    build_model(model_equations::AbstractString,R;df::Float64=4.0)

* Build models from **model equations** with residual varainces **R** and degree
of freedom for residual variance **df** defaults to 4.0.

```julia
#single-trait
model_equations = "BW = intercept + age + sex"
R               = 6.72
models          = build_model(model_equations,R);

#multi-trait
model_equations = "BW = intercept + age + sex;
                   CW = intercept + age + sex";
R               = [6.72   24.84
                   24.84  708.41]
models          = build_model(model_equations,R);
```
"""
function build_model(model_equations::AbstractString,R;df=4)
    if !(typeof(model_equations)<:AbstractString) || model_equations==""
        error("Model equations are wrong.\n
        An example for single-trait analyses is
        \"bw = age + herd\"\n
        An example for multi-trait analyses is
        \"bw = age + herd;
         cw = age + herd\" \n")
    end

    modelVec   = split(model_equations,[';','\n'],keep=false)
    nModels    = size(modelVec,1)
    lhsVec     = Symbol[] #:y
    modelTerms = ModelTerm[] #initiate outside for loop
    dict       = Dict{AbstractString,ModelTerm}()
    for (m,model) = enumerate(modelVec)
        lhsRhs = split(model,"=")                  #"y2","A+B+A*B"
        lhsVec = [lhsVec;symbol(strip(lhsRhs[1]))] #:y2
        rhsVec = split(strip(lhsRhs[2]),"+")       #"A","B","A*B"
        mTrms  = [ModelTerm(trmStr,m) for trmStr in rhsVec]
        modelTerms  = [modelTerms;mTrms] #vector of ModelTerm
        for (i,trm) = enumerate(modelTerms)
            dict[trm.trmStr] = modelTerms[i]
        end
    end
    return MME(nModels,modelVec,modelTerms,dict,lhsVec,map(Float64,R),Float64(df))
end


"""
    set_covariate(mme::MME,covStr::AbstractString...)

* set variables as covariates; covStr is a string of variables


```julia
#After running build_model, two ways to set the variable age and year as covariates are as follows.
set_covariate(models,"age","year")
set_covariate(models,"age year")
```
"""
function set_covariate(mme::MME,covStr::AbstractString...)
    covVec=[]
    for i in covStr
        covVec = [covVec;split(i," ",keep=false)]
    end
    mme.covVec = [mme.covVec;[Symbol(i) for i in covVec]]
    #@printf("Variables %-10s are set to be covariates.",mme.covVec)
end

#fill up str and val for each ModelTerm
function getData(trm::ModelTerm,df::DataFrame,mme::MME) #ModelTerm("1:A*B")
    nObs    = size(df,1)
    trm.str = Array(AbstractString,nObs)
    trm.val = Array(Float64,nObs)

    if trm.factors[1] == :intercept
        str = fill("intercept",nObs)
        val = fill(1.0,nObs)
    else #for ModelTerm object such as "1:A" or "1:A*B"
        myDf = df[trm.factors]                          #:A,:B
        if trm.factors[1] in mme.covVec #covariate
            str = fill(string(trm.factors[1]),nObs)     #["A","A",...]
            val = df[trm.factors[1]]                    #df[:A]
        else #factors (also animal or maternal effects)
            str = [string(i) for i in df[trm.factors[1]]] #["A1","A2","A1",...]
            val = fill(1.0,nObs)
        end

        #for ModelTerm object such as "A*B" whose nFactors>1
        for i=2:trm.nFactors
            if trm.factors[i] in mme.covVec
                str = str .* fill(" * "*string(trm.factors[i]),nObs) #["A * B","A * B",...]
                val = val .* df[trm.factors[i]]
            else
                str = str .* fill(" * ",nObs) .* [string(j) for j in df[trm.factors[i]]]
                val = val .* fill(1.0,nObs)  #["A1 * B1","A2 * B2",...]
            end
        end
    end

    trm.str = str
    trm.val = val
end

getFactor1(str) = [strip(i) for i in split(str,"*")][1] #using in may be better. maybe age*animal
                                                        #Bug: can only use animal*age, not age*animal
#make incidence matrix for each ModelTerm
function getX(trm::ModelTerm,mme::MME)
    pedSize = 0
    nObs  = size(trm.str,1)
    if trm.trmStr in mme.pedTrmVec #random variables(pedigree),e.g."Animal","Animal*age"
        trm.names   = PedModule.getIDs(mme.ped)
        trm.nLevels = length(mme.ped.idMap)
        xj          = round(Int64,[mme.ped.idMap[getFactor1(i)].seqID for i in trm.str[trm.str .!= "0"]])#remove founder(fit maternal effects)
    else
        dict,trm.names  = mkDict(trm.str)
        trm.nLevels     = length(dict)
        xj              = round(Int64,[dict[i] for i in trm.str]) #column index
    end
    xi    = (trm.iModel-1)*nObs + collect(1:nObs)  #row index
    xv    = trm.val                                #value

    #remove "0",introducing by maternal effects
    xv    = xv[trm.str .!= "0"]
    xi    = xi[trm.str .!= "0"]

    #some animal ID may be missing in data (df),
    #below to ensure number of columns for BV = number of animals
    if mme.ped!=0
        pedSize = length(mme.ped.idMap)
        if trm.trmStr in mme.pedTrmVec
            # This is to ensure the X matrix for
            # additive effect has the correct number of columns
            ii = 1         # adding a zero to
            jj = pedSize   # the last column in row 1
            vv = [0.0]
            xi = [xi;ii]
            xj = [xj;jj]
            xv = [xv;vv]
        end
    end
    #ensure X has nObs*nModels rows
    nModels = size(mme.lhsVec,1)
    xi = [xi;1;nObs*nModels] #if (1,1) and (nObs*nModels,1) already exist
    xj = [xj;1;1]            # add 0 to X
    xv = [xv;0;0]
    trm.X = sparse(xi,xj,xv)
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

    #make incidence matrices for each term
    for trm in mme.modelTerms
        getData(trm,df,mme)
        getX(trm,mme)
    end
    n   = size(mme.modelTerms,1)
    trm = mme.modelTerms[1]
    X   = trm.X
    #concat incidence matrix for each term
    for i=2:n
        trm = mme.modelTerms[i]
        X = [X trm.X]
    end

    #make response vector
    y = convert(Array,df[mme.lhsVec[1]],0.0) #NA to zero
    for i=2:size(mme.lhsVec,1)
        y    = [y; convert(Array,df[mme.lhsVec[i]],0.0)]
    end
    nInd  = size(y,1)
    ii    = 1:nInd
    jj    = fill(1,nInd)
    vv    = y
    ySparse = sparse(ii,jj,vv)

    #make lhs and rhs for mixed model equations
    mme.X       = X
    mme.ySparse = ySparse
    if mme.nModels>1 #multi-trait
      Ri         = mkRi(mme,df)
      mme.mmeLhs = X'Ri*X
      mme.mmeRhs = X'Ri*ySparse
    elseif mme.nModels==1 #single-trait (lambda version)
      mme.mmeLhs = X'X
      mme.mmeRhs = X'ySparse
    end

    if mme.ped != 0
        ii,jj,vv = PedModule.HAi(mme.ped)
        HAi = sparse(ii,jj,vv)
        mme.Ai = HAi'HAi
        addA(mme::MME)
    end

    #iid random effects,NEED another addlambda for multi-trait
    if mme.nModels==1 #single-trait
      addLambdas(mme)
    end
end
