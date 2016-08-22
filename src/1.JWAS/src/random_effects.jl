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

#set **specific ModelTerm** to random
"""
    set_random(mme::MME,randomStr::AbstractString,ped::Pedigree, G;df=4)

* set variables as random polygenic effects with pedigree information **ped**, variances **G** whose degree of freedom **df** defaults to 4.0.

```julia
model_equations = "BW = intercept + age + sex + Animal;
                   CW = intercept + age + sex + Animal";
model           = build_model(model_equations,R);
ped             = get_pedigree(pedfile);
G               = [6.72   2.84
                   2.84  8.41]
set_random(model,"Animal", ped,G)
```
"""
function set_random(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G;df=4)
    pedTrmVec = split(randomStr," ",keep=false)  # "animal" or "animal animal*age"

    #add model number => "1:animal"
    res = []
    for trm in pedTrmVec #"animal","animal animal*age"
        for (m,model) = enumerate(mme.modelVec)
            strVec  = split(model,['=','+'])
            strpVec = [strip(i) for i in strVec]
            if trm in strpVec
                res = [res;string(m)*":"*trm]
            end
        end
    end #"1:animal","1:animal animal*age"
    if length(res) != size(G,1)
        error("Dimensions must match. G should be a ",length(res)," x ",length(res)," matrix.\n")
    end
    mme.pedTrmVec = res
    mme.ped = ped

    if mme.nModels!=1 #multi-trait
      mme.Gi = inv(G)
    else #single-trait
      if issubtype(typeof(G),Number)==true #convert scalar G to 1x1 matrix
          G=reshape([G],1,1)
      end
      mme.GiOld = zeros(G)
      mme.GiNew = inv(G)
    end

    mme.df.polygenic=Float64(df)
    nothing
end

#May not proper, assume no covariance between traits for iid random #not require wishart distribution
"""
    set_random(mme::MME,randomStr::AbstractString,vc::Float64; df::Float64))

set variables as iid random effects
"""
function set_random(mme::MME,randomStr::AbstractString, vc::Float64; df=4)
    randTrmVec = split(randomStr," ",keep=false)  # "herd"
    res = []
    for trm in randTrmVec #add model number => "1:herd"
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
      randomEffect = RandomEffect(trm,1.0,vc,df,scale,zeros(1,1))#ROld/vcOld=0
      push!(mme.rndTrmVec,randomEffect)
    end

    mme.df.random=Float64(df)
    nothing
end

function addA(mme::MME)
    pedTrmVec = mme.pedTrmVec
    for (i,trmi) = enumerate(pedTrmVec)
        pedTrmi  = mme.modelTermDict[trmi]
        startPosi  = pedTrmi.startPos
        endPosi    = startPosi + pedTrmi.nLevels - 1
        for (j,trmj) = enumerate(pedTrmVec)
            pedTrmj  = mme.modelTermDict[trmj]
            startPosj= pedTrmj.startPos
            endPosj  = startPosj + pedTrmj.nLevels - 1
            #lamda version or not
            myaddA=(mme.nModels!=1)?(mme.Ai*mme.Gi[i,j]):(mme.Ai*(mme.GiNew[i,j]*mme.RNew - mme.GiOld[i,j]*mme.ROld))
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
        if mme.nModels==1
          lambdaDiff = mme.RNew/effect.vcNew - mme.ROld/effect.vcOld
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] =
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] + speye(trmi.nLevels)*lambdaDiff
        else
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] =
          mme.mmeLhs[startPosi:endPosi,startPosi:endPosi] + speye(trmi.nLevels)*(1/effect.vcNew)
        end
    end
end
