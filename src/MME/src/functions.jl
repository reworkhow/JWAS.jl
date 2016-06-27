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

function getNames(mme)
    names = Array(AbstractString,0)
    for trm in mme.modelTerms
        for name in trm.names
            push!(names,trm.trmStr*": "*name)
        end
    end
    return names
end

function getTerm(trmStr)
    trm = ModelTerm(trmStr,0,[],[],[],0,0,spzeros(0,0),[])
    factorVec = split(trmStr,"*")
    trm.nFactors = length(factorVec)
    trm.factors = [symbol(strip(f)) for f in factorVec]
    return trm
end

function initMME(modelEquation::AbstractString,R::Float64)
    if modelEquation==""
        error("modelEquation is empty\n")
    end
    lhsRhs = split(modelEquation,"=") #"y","A+B+A*B"
    lhs = symbol(strip(lhsRhs[1]))    #:y
    rhs = strip(lhsRhs[2])            #"A+B+A*B"
    rhsVec = split(rhs,"+")           #"A","B","A*B" 
    dict = Dict{AbstractString,ModelTerm}()
    modelTerms = [getTerm(strip(trmStr)) for trmStr in rhsVec]
    for (i,trm) = enumerate(modelTerms)
        dict[trm.trmStr] = modelTerms[i] #modelTermDict::Dict{AbstractString,ModelTerm}
    end
    return MME(modelEquation,modelTerms,dict,lhs,R)
end

function getData(trm::ModelTerm,df::DataFrame,mme::MME)  #ModelTerm("A*B")
    nObs = size(df,1)
    trm.str = Array(AbstractString,nObs)
    trm.val = Array(Float64,nObs)

    if trm.factors[1] == :intercept
        str = fill(string(trm.factors[1]),nObs)
        val = fill(1.0,nObs)
        trm.str = str
        trm.val = val
        return
    end

    myDf = df[trm.factors]   #:A,:B

    #for ModelTerm object such as "A" or "A*B"
    if trm.factors[1] in mme.covVec
        str = fill(string(trm.factors[1]),nObs)       #["A","A",...]
        val = df[trm.factors[1]]                      #df[:A]
    else
        str = [string(i) for i in df[trm.factors[1]]] #["A1","A2","A1",...]
        val = fill(1.0,nObs)
    end

    #for ModelTerm object such as "A*B" whose nFactors>1 
    for i=2:trm.nFactors
        if trm.factors[i] in mme.covVec
            str = str .* fill(" x "*string(trm.factors[i]),nObs) #["A x B","A x B",...]
            val = val .* df[trm.factors[i]] 
        else
            str = str .* fill(" x ",nObs) .* [string(j) for j in df[trm.factors[i]]] 
            val = val .* fill(1.0,nObs)                          #["A1 X B1","A2 x B2",...]                
        end
    end
    
    trm.str = str #detailed example in types.jl
    trm.val = val
end

function covList(mme::MME, covStr::AbstractString) #set variable as covariate
    covVec = split(covStr," ",keep=false)
    mme.covVec = [symbol(i) for i in covVec]
    nothing
end

getFactor1(str) = [strip(i) for i in split(str,"x")][1] #using in may be better. maybe age*animal

function getX(trm::ModelTerm,mme::MME) #make incidence matrix
    pedSize = 0
    nObs  = size(trm.str,1)
    if trm.trmStr in mme.pedTrmVec   
        #random variables related to A(pedigree),e.g."Animal"
        trm.names   = PedModule.getIDs(mme.ped)
        trm.nLevels = length(mme.ped.idMap)
        xj = round(Int64,[mme.ped.idMap[getFactor1(i)].seqID for i in trm.str]) #column index
    else   
        #fixed effects
        dict,trm.names  = mkDict(trm.str)
        trm.nLevels     = length(dict)
        xj    =  round(Int64,[dict[i] for i in trm.str]) #column index
    end
    xi    = 1:nObs  #row index
    xv    = trm.val #value
    
    #some animal ID may be missing in data (df), 
    #below to ensure number of columns for BV = number of animals 
    if mme.ped!=0
        pedSize = length(mme.ped.idMap)
        if trm.trmStr in mme.pedTrmVec
            # adding a zero to the last column in row 1
            ii = 1         
            jj = pedSize
            vv = [0.0]
            xi = [xi;ii]
            xj = [xj;jj]
            xv = [xv;vv]
        end
    end
    trm.X = sparse(xi,xj,xv)
    trm.startPos = mme.mmePos
    mme.mmePos  += trm.nLevels
end

function getMME(mme::MME, df::DataFrame)
    for trm in mme.modelTerms
        #make incidence matrices for each term
        getData(trm,df,mme)
        getX(trm,mme)
    end
    n   = size(mme.modelTerms,1)
    trm = mme.modelTerms[1]
    X   = trm.X
    for i=2:n
        #get incidence matrix for all
        trm = mme.modelTerms[i]
        X = [X trm.X]
    end
    y    = df[mme.lhs]
    nObs = size(y,1)
    ii = 1:nObs
    jj = fill(1,nObs)
    vv = y
    nRowsX = size(X,1)
    if nRowsX > nObs  ###?????? nRowsX=nObs #should delete
        ii = [ii,nRowsX]
        jj = [jj,1]
        vv = [vv,0.0]
    end
    ySparse = sparse(ii,jj,vv)
    mme.X = X
    mme.ySparse = ySparse
    mme.mmeLhs = X'X
    mme.mmeRhs = X'ySparse
    if mme.ped != 0
        ii,jj,vv = PedModule.HAi(mme.ped)#cholesky??
        HAi = sparse(ii,jj,vv)
        mme.Ai = HAi'HAi
        addA(mme::MME)
    end
    addLambdas(mme)
end

function setAsRandom(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G)
    mme.pedTrmVec = split(randomStr," ",keep=false)
    mme.ped = ped
    if issubtype(typeof(G),Number)==true  
        #when G is a scalar instead of matrix, make G 1x1 Matrix
        G=reshape([G],1,1)  
    end
    mme.GiOld = zeros(G)
    mme.GiNew = inv(G)
    nothing
end

function setAsRandom(mme::MME,randomStr::AbstractString, vc::Float64, df::Float64)
    trm  = mme.modelTermDict[randomStr]
    scale = vc*(df-2)/df
    randomEffect = RandomEffect(trm,1.0,vc,df,scale,Array(Float64,1))#ROld/vcOld=0
    push!(mme.rndTrmVec,randomEffect)
end

function addA(mme::MME) #add Ainv*lambda
    pedTrmVec = mme.pedTrmVec
    for (i,trmi) = enumerate(pedTrmVec)
        pedTrmi  = mme.modelTermDict[trmi]
        startPosi  = pedTrmi.startPos
        endPosi    = startPosi + pedTrmi.nLevels - 1
        for (j,trmj) = enumerate(pedTrmVec)
            pedTrmj  = mme.modelTermDict[trmj]
            startPosj  = pedTrmj.startPos
            endPosj    = startPosj + pedTrmj.nLevels - 1
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] =
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] +             
            mme.Ai*(mme.GiNew[i,j]*mme.RNew - mme.GiOld[i,j]*mme.ROld)
        end
    end
end

function addLambdas(mme::MME) #for iid random effects
    for effect in  mme.rndTrmVec
        trmi       = effect.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        lambdaDiff = mme.RNew/effect.vcNew - mme.ROld/effect.vcOld
        mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] = 
        mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] + speye(trmi.nLevels)*lambdaDiff
    end   
end

function getSolJ(mme::MME, df::DataFrame)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) Jacobi(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,0.3,tol=0.000001)]
end

function getSolG(mme::MME, df::DataFrame;outFreq=10)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) GaussSeidel(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,tol=0.000001,output=outFreq)]
end

function getSolGibbs(mme::MME, df::DataFrame;nIter=50000,outFreq=100)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    return [getNames(mme) Gibbs(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,mme.RNew,nIter,outFreq=outFreq)]
end
