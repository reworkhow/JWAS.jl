"""
    get_pedigree(pedfile::AbstractString)
* Get pedigree informtion from a pedigree file.
* File format:

```
a 0 0
b 0 0
c a b
d a c
```
"""
function get_pedigree(pedfile::AbstractString)
  PedModule.mkPed(pedfile)
end

################################################################################
#set specific ModelTerm to random
#e.g. Animal, Animal*Age, Maternal
################################################################################
"""
    set_random(mme::MME,randomStr::AbstractString,ped::Pedigree, G;df=4)

* set variables as random polygenic effects with pedigree information **ped**, variances **G** whose degree of freedom **df** defaults to 4.0.

```julia
#single-trait (example 1)
model_equation  = "y = intercept + Age + Animal"
model           = build_model(model_equation,R)
ped             = get_pedigree(pedfile)
G               = 1.6
set_random(model,"Animal Animal*Age", ped,G)

#single-trait (example 2)
model_equation  = "y = intercept + Age + Animal + Animal*Age"
model           = build_model(model_equation,R)
ped             = get_pedigree(pedfile)
G               = [1.6   0.2
                   0.2  1.0]
set_random(model,"Animal Animal*Age", ped,G)

#multi-trait
model_equations = "BW = intercept + age + sex + Animal
                   CW = intercept + age + sex + Animal"
model           = build_model(model_equations,R);
ped             = get_pedigree(pedfile);
G               = [6.72   2.84
                   2.84  8.41]
set_random(model,"Animal", ped,G)
```
"""
function set_random(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G;df=4)
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
          error(trm," is not found in model equations.")
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

"""
    set_random(mme::MME,randomStr::AbstractString,G;df=4)

* set variables as i.i.d random effects with variances **G** whose degree of freedom **df** defaults to 4.0.

```julia
#single-trait (example 1)
model_equation  = "y = intercept + litter + sex"
model           = build_model(model_equation,R)
G               = 0.6
set_random(model,"litter",G)

#multi-trait
model_equations = "BW = intercept + litter + sex
                   CW = intercept + litter + sex"
model           = build_model(model_equations,R);
G               = [3.72  1.84
                   1.84  3.41]
set_random(model,"litter",G)
```
"""
function set_random(mme::MME,randomStr::AbstractString, G; df=4)
    G = map(Float64,G)
    randTrmVec = split(randomStr," ",keep=false)  # "herd"
    res = []
    for trm in randTrmVec #add model number => "1:herd"
      for (m,model) = enumerate(mme.modelVec)
        strVec  = split(model,['=','+'])
        strpVec = [strip(i) for i in strVec]
        if trm in strpVec
          res = [res;string(m)*":"*trm]
        else
          error(trm," is not found in model equations.")
        end
      end
    end #"1:herd","2:herd"
    if length(res) != size(G,1)
      error("Dimensions must match. The covariance matrix (G) should be a ",length(res)," x ",length(res)," matrix.\n")
    end

    if mme.nModels!=1 #multi-trait
      mme.Gi = inv(G)
    else              #single-trait
      if issubtype(typeof(G),Number)==true #convert scalar G to 1x1 matrix
        G=reshape([G],1,1)
      end
      mme.GiOld = zeros(G)
      mme.GiNew = inv(G)
    end


    for term_string in res
      trm  = mme.modelTermDict[term_string]
      scale = vc*(df-2)/df
      randomEffect = RandomEffect(trm,1.0,vc,df,scale,zeros(1,1))#ROld/vcOld=0
      push!(mme.rndTrmVec,randomEffect)
    end

    mme.df.random=Float64(df)
    nothing
end
################################################################################
#*******************************************************************************
#following facts that scalar Inverse-Wishart(ν,S) = Inverse-Gamma(ν/2,S/2)=    *
#scale-inv-chi2(ν,S/ν), variances for random effects(non-marker) will be be    *
#sampled from Inverse-Wishart for coding simplicity thus prior for scalar      *
#variance is treated as a 1x1 matrix.                                          *
#*******************************************************************************
################################################################################


#construct MME for pedigree-based random effects part
function addA(mme::MME) #two terms,e.g.,"animal" and "maternal" may not near in MME
    pedTrmVec = mme.pedTrmVec
    for (i,trmi) = enumerate(pedTrmVec)
        pedTrmi  = mme.modelTermDict[trmi]
        startPosi  = pedTrmi.startPos
        endPosi    = startPosi + pedTrmi.nLevels - 1
        for (j,trmj) in enumerate(pedTrmVec)
            pedTrmj  = mme.modelTermDict[trmj]
            startPosj= pedTrmj.startPos
            endPosj  = startPosj + pedTrmj.nLevels - 1
            #lamda version (single trait) or not (multi-trait)
            myaddA   = (mme.nModels!=1)?(mme.Ai*mme.Gi[i,j]):(mme.Ai*(mme.GiNew[i,j]*mme.RNew - mme.GiOld[i,j]*mme.ROld))
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] =
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] + myaddA
        end
    end
end

################################################################################
#construct MME for iid random effects part
#single-trait: lamda version, assume effects are independent, e.g., litter and groups
#multi-trait: same effects in 2 traits are correlated,e.g., 1:litter, 2:litter
################################################################################
function addLambdas(mme::MME)
    for effect in  mme.rndTrmVec
        trmi       = effect.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        if mme.nModels==1
          lambdaDiff = mme.RNew/effect.vcNew - mme.ROld/effect.vcOld
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] =
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] + speye(trmi.nLevels)*lambdaDiff
        #else
        #  mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] =
        #    mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] + speye(trmi.nLevels)*(1/effect.vcNew)
        end
    end
end
