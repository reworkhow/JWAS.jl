################################################################################
#Get names for variables in order
################################################################################
function getNames(mme::MME)
    names = Array(AbstractString,0)
    for trm in mme.modelTerms
        for name in trm.names
            push!(names,trm.trmStr*" : "*name)
        end
    end
    return names
end

function getNames(trm::ModelTerm)
    names = Array(AbstractString,0)
    for name in trm.names
        push!(names,trm.trmStr*" : "*name)
    end
    return names
end

################################################################################
#Print out MME (no markers)
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
