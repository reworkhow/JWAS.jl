#the column name for group records has to be "group_records"
function get_T(data)
    uID = data[:,"group_records"] # group info
    yID = unique(uID) #unique groups

    #get an incidence matrix Z to reorder uID to yID by yID = Z*uID
    T = mkmat_incidence_factor(uID,yID)'
    return T
end

function group_records_setup(mme,df)
    if mme.nModels != 1
        error("Group records analysis is only supported in single-trait analysis!")
    end
    T = map(Float32,get_T(df))

    #update Incidence Matrix
    for i = 1:length(mme.modelTerms) #modify incidence matrix for non-genotype effects
        mme.modelTerms[i].X=T*mme.modelTerms[i].X
    end
    #concatenate all terms
    X   = mme.modelTerms[1].X
    for i=2:length(mme.modelTerms)
       X = [X mme.modelTerms[i].X]
    end
    mme.X       = X

    #update phenotypes
    mme.ySparse = T*mme.ySparse #modify phenotypes

    #heterogeneous residuals
    invweights     = 1 ./ convert(Array,diag(T*T'))
    mme.invweights = (mme.MCMCinfo == false || mme.MCMCinfo.double_precision ? Float64.(invweights) : Float32.(invweights))

    #update Mixed Model Equation
    mme.mmeLhs = X'*Diagonal(mme.invweights)*X
    mme.mmeRhs = X'*Diagonal(mme.invweights)*mme.ySparse
    #Random effects parts in MME
    #random_term.GiNew*mme.R - random_term.GiOld*mme.ROld
    for random_term in mme.rndTrmVec #trick
      random_term.GiOld = zero(random_term.GiOld)
    end
    addVinv(mme)
    for random_term in mme.rndTrmVec #trick
      random_term.GiOld = copy(random_term.GiNew)
    end
    dropzeros!(mme.mmeLhs)
    dropzeros!(mme.mmeRhs)

    #update genotypes
    if mme.M != 0
        for Mi in mme.M #modify incidence matrix for genotype effects
            Mi.genotypes = T*Mi.genotypes
            #Mi.obsID     = mme.obsID
            #Mi.nObs      = length(mme.obsID)
        end
    end
end
