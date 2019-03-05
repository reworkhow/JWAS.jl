################################################################################
#Set specific ModelTerm as random using pedigree information
#e.g. Animal, Animal*Age, Maternal
################################################################################
"""
    set_random(mme::MME,randomStr::AbstractString,ped::Pedigree, G;df=4)

* set variables as random polygenic effects with pedigree information **ped** and
  variances **G** with degree of freedom **df**, defaulting to `4.0`.

```julia
#single-trait (example 1)
model_equation  = "y = intercept + age + animal"
model           = build_model(model_equation,R)
ped             = get_pedigree(pedfile)
G               = 1.6
set_random(model,"animal", ped, G)

#single-trait (example 2)
model_equation  = "y = intercept + age + animal + animal*age"
model           = build_model(model_equation,R)
ped             = get_pedigree(pedfile)
G               = [1.6   0.2
                   0.2  1.0]
set_random(model,"animal animal*age", ped,G)

#multi-trait
model_equations = "BW = intercept + age + sex + animal
                   CW = intercept + age + sex + animal"
model           = build_model(model_equations,R);
ped             = get_pedigree(pedfile);
G               = [6.72   2.84
                   2.84  8.41]
set_random(model,"animal",ped,G)
```
"""
function set_random(mme::MME,randomStr::AbstractString,ped::PedModule.Pedigree, G;df=4)
    if !isposdef(G)
        error("The covariance matrix is not positive definite.")
    end

    pedTrmVec = split(randomStr," ",keepempty=false)  # "animal animal*age"

    #add model equation number to variables;
    #"animal" => "1:animal"; "animal*age"=>"1:animal*age"
    res = []
    for trm in pedTrmVec
      for (m,model) = enumerate(mme.modelVec)
        strVec  = split(model,['=','+'])
        strpVec = [strip(i) for i in strVec]
        if trm in strpVec
          res = [res;string(m)*":"*trm]
        else
          printstyled(trm," is not found in model equation ",string(m),".\n",bold=false,color=:red)
        end
      end
    end
    if length(res) != size(G,1)
      error("Dimensions must match. The covariance matrix (G) should be a ",length(res)," x ",length(res)," matrix.\n")
    end
    mme.pedTrmVec = res
    if mme.ped == 0
        mme.ped = deepcopy(ped)
    else
        error("Pedigree information is already added.")
    end


    if typeof(G)<:Number #single variable single-trait, convert scalar G to 1x1 matrix
      G=reshape([G],1,1)
    end
    mme.Gi    = inv(G) #multi-trait
    mme.GiOld = inv(G) #single-trait (maybe multiple correlated genetic effects
    mme.GiNew = inv(G) #single-trait (in single-trait

    mme.df.polygenic=Float64(df)

    nothing
end
############################################################################
#Set specific ModelTerm as random  (Not specific for IID or A(pedigree))
#useful for imputaion residual in single-step methods (var(ϵ)^{-1}=A^{nn} )
############################################################################
"""
    set_random(mme::MME,randomStr::AbstractString,G;Vinv=0,names=[],df=4)

* set variables as random effects, defaulting to i.i.d effects, with variances **G** with
  degree of freedom **df**, defaulting to 4.0.
* the random effects are assumed to be i.i.d by default and it can be defined with any
  (inverse of) covariance structure **Vinv** with its index (row names) provided by **names**.

```julia
#single-trait (i.i.d randome effects)
model_equation  = "y = intercept + litter + sex"
model           = build_model(model_equation,R)
G               = 0.6
set_random(model,"litter",G)

#multi-trait (i.i.d randome effects)
model_equations = "BW = intercept + litter + sex
                   CW = intercept + litter + sex"
model           = build_model(model_equations,R);
G               = [3.72  1.84
                   1.84  3.41]
set_random(model,"litter",G)

#single-trait (randome effects with specific covariance structures)
model_equation  = "y = intercept + litter + sex"
model           = build_model(model_equation,R)
V               = [1.0  0.5 0.25
                   0.5  1.0 0.5
                   0.25 0.5 1.0]
G               = 0.6
set_random(model,"litter",G,Vinv=inv(V),names=[a1;a2;a3])
```
"""
function set_random(mme::MME,randomStr::AbstractString,G;Vinv=0,names=[],df=4.0)
    if !isposdef(G)
        error("The covariance matrix is not positive definite.")
    end
    G = map(Float64,G)
    df= Float64(df)
    randTrmVec = split(randomStr," ",keepempty=false)  # "litter"

    #add model equation number to variables; "litter"=>"1:litter";"ϵ"=>"1:ϵ"
    for trm in randTrmVec
      res = []
      for (m,model) = enumerate(mme.modelVec)
        strVec  = split(model,['=','+'])
        strpVec = [strip(i) for i in strVec]
        if trm in strpVec || trm == "ϵ"
          mtrm= string(m)*":"*trm
          res = [res;mtrm]
          #*********************
          #phenotype IDs is a subset of names (Vinv)
          mme.modelTermDict[mtrm].names=names
          #*********************
        else
          printstyled(trm," is not found in model equation ",string(m),".\n",bold=false,color=:red)
        end
      end
      if length(res) != size(G,1)
        error("Dimensions must match. The covariance matrix (G) should be a ",length(res)," x ",length(res)," matrix.\n")
      end
      if typeof(G)<:Number ##single variable single-trait, convert scalar G to 1x1 matrix
        G=reshape([G],1,1)
      end
      Gi    = inv(G)    #multi-trait
      GiOld = inv(G)    #single-trait (maybe multiple correlated random effects
      GiNew = inv(G)    #single-trait (thus G is a matrix

      term_array   = res
      df           = df+length(term_array)
      scale        = G*(df-length(term_array)-1)  #G*(df-2)/df #from inv χ to inv-wishat
      randomEffect = RandomEffect(term_array,Gi,GiOld,GiNew,df,scale,Vinv,names)
      push!(mme.rndTrmVec,randomEffect)
    end
    nothing
end



################################################################################
#ADD TO MIXED MODEL EQUATIONS FOR THE RANDOM EFFECTS PARTS
################################################################################
#NOTE: SINGLE TRAIT
#The equation Ai*(GiNew*RNew - GiOld*ROld) is used to update Ai part in LHS
#The 1st time to add Ai to set up MME,
#mme.GiOld == zeros(G),mme.GiNew == inv(G),mme.Rnew == mme.Rold= R
#After that, if variances are constant,
#mme.GiOld == mme.GiNew; mme.Rnew == mme.Rold
#If sampling genetic variances, mme.Ginew is updated with a new sample, then
#LHS is update as Ai*(GiNew*RNew - GiOld*ROld), then GiOld = Ginew
#if sample residual variances, similar approaches to update the LHS is used.
################################################################################

#Construct MME for pedigree-based random effects part
function addA(mme::MME) #two terms,e.g.,"animal" and "maternal" may not be near in MME
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
            myaddA   = (mme.nModels!=1) ? (mme.Ai*mme.Gi[i,j]) : (mme.Ai*(mme.GiNew[i,j]*mme.RNew - mme.GiOld[i,j]*mme.ROld))
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] =
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] + myaddA
        end
    end
end

################################################################################
#Construct MME for random effects with specific covariance structures (including iid )
#Single-trait:
#lamda version, e.g., litter and groups
#1) use `set_random(model,"1:A 1:B",G)` to estimate cov(1:A,1:B) unless A.names == B.names,
#2) use `set_random(model,"1:A",G1);`set_random(model,"1:B",G2)` to make 1:A and 1:B independent.
#
#Multi-trait:
#e.g., 1:litter, 2:litter
################################################################################
function addLambdas(mme::MME)
    for random_term in mme.rndTrmVec
      term_array = random_term.term_array
      Vi         = (random_term.Vinv!=0) ? random_term.Vinv : SparseMatrixCSC{Float64}(I, mme.modelTermDict[term_array[1]].nLevels, mme.modelTermDict[term_array[1]].nLevels)
      for (i,termi) = enumerate(term_array)
          randTrmi   = mme.modelTermDict[termi]
          startPosi  = randTrmi.startPos
          endPosi    = startPosi + randTrmi.nLevels - 1
          for (j,termj) in enumerate(term_array)
            randTrmj    = mme.modelTermDict[termj]
            startPosj   = randTrmj.startPos
            endPosj     = startPosj + randTrmj.nLevels - 1
            myaddLambda = (mme.nModels!=1) ? (Vi*random_term.Gi[i,j]) : (Vi*(random_term.GiNew[i,j]*mme.RNew - random_term.GiOld[i,j]*mme.ROld))
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] =
            mme.mmeLhs[startPosi:endPosi,startPosj:endPosj] + myaddLambda
          end
      end
    end
end
