function censored_trait_setup!(mme,lower_bound,upper_bound)
    nInd                 = length(mme.obsID)
    nTrait               = mme.nModels
    censored_trait_index = mme.censored_trait_index #e.g., [2,4] means the 2nd and 4th traits are censored
    n_censored_trait     = length(censored_trait_index)


    ySparse=reshape(Matrix(mme.ySparse),nInd,nTrait)

    lower_bound[:,censored_trait_index]  = Float64.(ySparse)[:,censored_trait_index]  # nInd-by-#censored_trait
    upper_bound[:,censored_trait_index]  = Float64.(mme.MCMCinfo.censored_trait)      # nInd-by-#censored_trait

    ############################################################################
    # liability (mme.ySparse)
    ############################################################################
    for t in 1:nTrait
        if mme.lhsTag[t]=="censored"
            for i in 1:nInd
                if lower_bound[i,t] == -Inf
                    ySparse[i,t] = upper_bound[i,t]
                elseif upper_bound[i,t] == Inf
                    ySparse[i,t] = lower_bound[i,t]
                else
                    ySparse[i,t] = rand(lower_bound[i,t]:upper_bound[i,t])
                end
            end
        end
    end
    mme.ySparse = reshape(ySparse,nInd*nTrait,1)
    return lower_bound,upper_bound
end



#e.g., censored trait is ["y1"].
# old: "y1   = intercept + genotypes; y11 = intercept + genotypes"
# new: "y1_l = intercept + genotypes; y11 = intercept + genotypes"
function censored_trait_model_equation(model_equations,censored_trait)
    for t in censored_trait
        model_equations = replace(model_equations, t*r"\s*=" => t*"_l=") # "\s*" means zero or more space;
    end #here the censored_trait must be in the left of "=" to keep "y11" unchanged.
    return model_equations
end
