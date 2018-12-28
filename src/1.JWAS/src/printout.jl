################################################################################
#Get names for variables in order
################################################################################
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

################################################################################
#Print out model information
################################################################################
"""
    getinfo(model::MME)

* Print out model information.
"""

#more details later
function getinfo(model;data=false)
  println("A Linear Mixed Model was build using model equations:\n")
  for i in model.modelVec
    println(i)
  end
  println()
  println("Model Information:\n")
  @printf("%-15s %-12s %-10s %11s\n","Term","C/F","F/R","nLevels")

  random_effects=Array{AbstractString,1}()
  if model.pedTrmVec != 0
    for i in model.pedTrmVec
        push!(random_effects,split(i,':')[end])
    end
  end
  for i in model.rndTrmVec
      for j in i.term_array
          push!(random_effects,split(j,':')[end])
      end
  end

  terms=[]
  for i in model.modelTerms
    term    = split(i.trmStr,':')[end]
    if term in terms
        continue
    else
        push!(terms,term)
    end

    nLevels = i.nLevels
    fixed   = (term in random_effects) ? "random" : "fixed"
    factor  = (nLevels==1) ? "covariate" : "factor"

    if size(model.mmeRhs)==()&&data==false
      nLevels="NA"
    elseif size(model.mmeRhs)==()
      getMME(model,data)
    end

    if term =="intercept"
        factor="factor"
    elseif length(split(term,'*'))!=1
        factor="interaction"
    end

    @printf("%-15s %-12s %-10s %11s\n",term,factor,fixed,nLevels)
  end
  println()
  #incidence matrix , #elements non-zero elements
  #"incomplete or complete",genomic data
end
