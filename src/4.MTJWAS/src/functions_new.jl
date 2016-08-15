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

function getTerm(trmStr,m)
    trm = ModelTerm(m,string(m)*":"*trmStr)
    factorVec = split(trmStr,"*")
    trm.nFactors = length(factorVec)
    trm.factors = [symbol(strip(f)) for f in factorVec]
    return trm
end

function initMME(models::AbstractString,R)
    # returns an MME object from the input string
    if models==""
        error("modelEquation is empty\n")
    end
    modelVec   = split(models,[';','\n'],keep=false)
    nModels    = size(modelVec,1)                    #:y
    lhsVec     = Symbol[]
    dict       = Dict{AbstractString,ModelTerm}()
    for (m,model) = enumerate(modelVec)
        lhsRhs = split(model,"=")                  #"y2","A+B+A*B"
        lhsVec = [lhsVec;symbol(strip(lhsRhs[1]))] #:y2
        rhsVec = split(strip(lhsRhs[2]),"+")       #"A","B","A*B"
        modelTerms  = [getTerm(strip(trmStr),m) for trmStr in rhsVec] #vector of ModelTerm
        for (i,trm) = enumerate(modelTerms)
            dict[trm.trmStr] = modelTerms[i]
        end
    end
    return MME(modelVec,modelTerms,dict,lhsVec,R)
end

function covList(mme::MME, covStr::AbstractString) #Better with multi-arguments
    covVec = split(covStr," ",keep=false)
    mme.covVec = [symbol(i) for i in covVec]
    nothing
end

#fill up str and val for each ModelTerm
function getData(trm::ModelTerm,df::DataFrame,mme::MME) #ModelTerm("A*B")
    nObs    = size(df,1)
    trm.str = Array(AbstractString,nObs)
    trm.val = Array(Float64,nObs)

    if trm.factors[1] == :intercept
        str = fill("intercept",nObs)
        val = fill(1.0,nObs)
    else #for ModelTerm object such as "A" or "A*B"
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
                str = str .* fill(" * "*string(trm.factors[i]),nObs) #["A x B","A x B",...]
                val = val .* df[trm.factors[i]]
            else
                str = str .* fill(" * ",nObs) .* [string(j) for j in df[trm.factors[i]]]
                val = val .* fill(1.0,nObs)  #["A1 X B1","A2 x B2",...]
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
    if trm.trmStr in mme.pedTrmVec #random variables related to A(pedigree),e.g."Animal"
        trm.names   = PedModule.getIDs(mme.ped)
        trm.nLevels = length(mme.ped.idMap)
        xj          = round(Int64,[mme.ped.idMap[getFactor1(i)].seqID for i in trm.str[trm.str .!= "0"]])#remove founder
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
    xj = [xj;1;1]            # add to it,(0)
    xv = [xv;0;0]
    trm.X = sparse(xi,xj,xv)
    trm.startPos = mme.mmePos
    mme.mmePos  += trm.nLevels
end

function getMME(mme::MME, df::DataFrame)
    #make incidence matrices for each term
    for trm in mme.modelTerms
        getData(trm,df,mme)
        getX(trm,mme)
    end
    n   = size(mme.modelTerms,1)
    trm = mme.modelTerms[1]
    X   = trm.X
    for i=2:n #concat incidence matrix for each term
        trm = mme.modelTerms[i]
        X = [X trm.X]
    end

    y = convert(Array,df[mme.lhsVec[1]],0.0) #NA to zero
    for i=2:size(mme.lhsVec,1)
        y    = [y; convert(Array,df[mme.lhsVec[i]],0.0)]
    end
    N  = size(y,1)
    ii = 1:N
    jj = fill(1,N)
    vv = y
    ySparse = sparse(ii,jj,vv)

    nObs = size(df,1)
    Ri = mkRi(mme,df)
    mme.X = X
    mme.ySparse = ySparse
    if length(mme.modelVec)!=1 #multi-trait
      mme.mmeLhs = X'Ri*X
      mme.mmeRhs = X'Ri*ySparse
    else #single-trait (lambda version)
      mme.mmeLhs = X'X
      mme.mmeRhs = X'ySparse
    end

    if mme.ped != 0
        ii,jj,vv = PedModule.HAi(mme.ped)
        HAi = sparse(ii,jj,vv)
        mme.Ai = HAi'HAi
        addA(mme::MME)
    end

    #NEED another addlambda for multi-trait
    if length(mme.modelVec)==1 #single-trait
      addLambdas(mme)
    end
end

#set **specific ModelTerm** to random
function setAsRandom(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G)
    pedTrmVec = split(randomStr," ",keep=false)  # "animal" or "animal animal*age"
    res = []
    for trm in pedTrmVec #add model number => "1:animal"
        for (m,model) = enumerate(mme.modelVec)
            strVec  = split(model,['=','+'])
            strpVec = [strip(i) for i in strVec]
            if trm in strpVec
                res = [res;string(m)*":"*trm]
            end
        end
    end
    mme.pedTrmVec = res
    mme.ped = ped

    if length(mme.modelVec)!=1 #multi-trait
      mme.Gi = inv(G)
    else #single-trait
      if issubtype(typeof(G),Number)==true #convert scalar G to 1x1 matrix
          G=reshape([G],1,1)
      end
      mme.GiOld = zeros(G)
      mme.GiNew = inv(G)
    end
    nothing
end

#May not proper, assume no covariance between traits for iid random #not require wishart distribution
function setAsRandom(mme::MME,randomStr::AbstractString, vc::Float64, df::Float64)
    randTrmVec = split(randomStr," ",keep=false)  # "herd"
    res = []
    for trm in randTrmVec #add model number => "1:animal"
        for (m,model) = enumerate(mme.modelVec)
            strVec  = split(model,['=','+'])
            strpVec = [strip(i) for i in strVec]
            if trm in strpVec
                res = [res;string(m)*":"*trm]
            end
        end
    end #"1:herd","2:herd"

    for term_string in res
      trm  = mme.modelTermDict[term_string]
      scale = vc*(df-2)/df
      randomEffect = RandomEffect(trm,1.0,vc,df,scale,Array(Float64,1))#ROld/vcOld=0
      push!(mme.rndTrmVec,randomEffect)
    end

    nothing
end

function setAsRandom(mme::MME,randomStr::AbstractString, vc::Float64, df::Float64)
    trm  = mme.modelTermDict[randomStr]
    scale = vc*(df-2)/df
    randomEffect = RandomEffect(trm,1.0,vc,df,scale,Array(Float64,1))#ROld/vcOld=0
    push!(mme.rndTrmVec,randomEffect)
end


function addA(mme::MME) #little different from single trait
    pedTrmVec = mme.pedTrmVec
    for (i,trmi) = enumerate(pedTrmVec)
        pedTrmi  = mme.modelTermDict[trmi]
        startPosi  = pedTrmi.startPos
        endPosi    = startPosi + pedTrmi.nLevels - 1
        for (j,trmj) = enumerate(pedTrmVec)
            pedTrmj  = mme.modelTermDict[trmj]
            startPosj  = pedTrmj.startPos
            endPosj    = startPosj + pedTrmj.nLevels - 1
            #lamda version or not
            myaddA=(length(mme.modelVec)!=1)?(mme.Ai*mme.Gi[i,j]):(mme.Ai*(mme.GiNew[i,j]*mme.RNew - mme.GiOld[i,j]*mme.ROld))
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] =
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] + myaddA
        end
    end
end

function addLambdas(mme::MME) #for iid random effects
    for effect in  mme.rndTrmVec
        trmi       = effect.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        if length(mme.modelVec)==1
          lambdaDiff = mme.RNew/effect.vcNew - mme.ROld/effect.vcOld
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] =
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] + speye(trmi.nLevels)*lambdaDiff
        else
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] =
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] + speye(trmi.nLevels)*(1/effect.vcNew)
        end
    end
end

function getNames(mme)
    names = Array(AbstractString,0)
    for trm in mme.modelTerms
        for name in trm.names
            push!(names,trm.trmStr*" : "*name)
        end
    end
    return names
end
