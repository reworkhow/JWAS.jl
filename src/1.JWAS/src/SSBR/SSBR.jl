type SSBR
    Ai_nn =Ai[1:num.pedn,1:num.pedn]
    Ai_ng = Ai[1:num.pedn,(num.pedn+1):num.ped]

end


function add_data_imputation_residual(mme,df)
#add a fake column for imputaion RESIDUAL
    IDs       = convert(Array,df[:,1])
    isnongeno = [ID in mme.ped.setNG for ID in IDs] #true/false
    df[:imputation_residual]=copy(df[:,1])
    df[:imputation_residual][!isnongeno]="0"
end

#modify model equations
"""
add to model an extra term: imputation_residual
"""
function add_term_imputation_residual(mme,G;df=4)
    for m in 1:nModels
        push!(mme.modelTerms,ModelTerm("imputation_residual",m))
        mme.dict[ModelTerm("imputation_residual",m).trmStr]=ModelTerm("imputation_residual",m)
    end


  modelVec   = [strip(i) for i in split(model_equations,[';','\n'],keep=false)]
  nModels    = size(modelVec,1)
  lhsVec     = Symbol[]    #:y, phenotypes
  modelTerms = ModelTerm[] #initialization outside for loop
  dict       = Dict{AbstractString,ModelTerm}()
  for (m,model) = enumerate(modelVec)
    lhsRhs = split(model,"=")                  #"y2","A+B+A*B"
    lhsVec = [lhsVec;Symbol(strip(lhsRhs[1]))] #:y2
    rhsVec = split(strip(lhsRhs[2]),"+")       #"A","B","A*B"
    mTrms  = [ModelTerm(strip(trmStr),m) for trmStr in rhsVec]
    modelTerms  = [modelTerms;mTrms]           #vector of ModelTerm
  end
  for (i,trm) = enumerate(modelTerms)          #make a dict for model terms
    dict[trm.trmStr] = modelTerms[i]
  end
  return MME(nModels,modelVec,modelTerms,dict,lhsVec,map(Float64,R),Float64(df))
  #ARE MODELVEC LHSVEC USED?
end



function calc_Ai(ped::PedModule.Pedigree,geno::Genotypes)
    num_pedn    = PedModule.genoSet!(geno.obsID,ped) #set order in A as geno id then non-geno id
    if mme.Ai == 0
        mme.Ai          = PedModule.AInverse(ped)
    else
        error("Something is wrong in SSBR.")
    end
    Ai_nn       = mme.Ai[1:num_pedn,1:num_pedn]
    Ai_ng       = mme.Ai[1:num_pedn,(num_pedn+1):size(mme.Ai,1)]
    return Ai_nn,Ai_ng
end

function make_MMats(geno::Genotypes,a_mats::AiMats,ped::PedModule.Pedigree)
    Mg   = Array{Float64,1}(geno.nObs,geno.nMarkers)
    MgID = Array{String,1}(geno.nObs)
    #reorder genotypes to get Mg with same order as Ai_gg
    for i in 1:geno.nObs
      id        = geno.obsID[i]
      row       = ped.idMap[id].seqID - num.pedn
      Mg[row,:] = geno.genotypes[i,:]
      MgID[row] = id
    end
    #Not store Mn
    #Mn = Ai_nn\(-Ai_ng*Mg)
    #M  = [Mn;Mg];
    mme.M.genotypes = Mg
    mme.M.obsID     = Mg

    Mfull = [Ai_nn\(-Ai_ng*Mg);Mg];
    IDs   = PedModule.getIDs(ped) #mme.M.obsID at the end
    gc()
end

function make_XWMats(jvecs,zmats,mmats,num::Numbers)#now fixed effects: μ
    Xn  = hcat(ones(num.yn), zmats.n*jvecs.n)
    Xg  = hcat(ones(num.yg), zmats.g*jvecs.g)
    X   =[Xn;
          Xg]

    #Wn = zmats.n*mmats.n
    #Wg = zmats.g*mmats.g
    #W  = [Wn;Wg];

    W  = zmats.full*mmats.full
    Wn = Array{Float64,2}(0,0)
    Wg = Array{Float64,2}(0,0)

    return XMats(X,Xn,Xg),WMats(W,Wn,Wg)
end

"""
    set_random(mme::MME,randomStr::AbstractString,Vinv,G;df=4)

* set variables as random polygenic effects with inverse of covariance matrix **Vinv** **ped**, variances **G** whose degree of freedom **df** defaults to 4.0.
"""
#For any random effects (Not specific for IID or A)
#useful for imputaion residual in single-step methods (var(ϵ)^{-1}=A^{nn} )
function set_random(mme::MME,randomStr::AbstractString,Vinv, G;df=4)
    pedTrmVec = split(randomStr," ",keep=false)  # "animal animal*age"

    #add model equation number; "animal" => "1:animal"
    res = []
    for trm in pedTrmVec
      for (m,model) = enumerate(mme.modelVec)
        strVec  = split(model,['=','+'])
        strpVec = [strip(i) for i in strVec]
        if trm in strpVec
          res = [res;string(m)*":"*trm]
        else
          info(trm," is not found in model equation ",string(m),".")
        end
      end
    end #"1:animal","1:animal*age"
    if length(res) != size(G,1)
      error("Dimensions must match. The covariance matrix (G) should be a ",length(res)," x ",length(res)," matrix.\n")
    end
    mme.pedTrmVec = res
    mme.ped = ped

    if mme.nModels!=1 #multi-trait
      mme.Gi = inv(G)
    else              #single-trait
      if issubtype(typeof(G),Number)==true #convert scalar G to 1x1 matrix
        G=reshape([G],1,1)
      end
      mme.GiOld = zeros(G)
      mme.GiNew = inv(G)
    end
    mme.df.polygenic=Float64(df)
    nothing
end
