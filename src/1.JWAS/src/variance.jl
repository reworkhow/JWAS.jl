#sample vairance of marker effects or residual
function sample_variance(x, n, df, scale)
    return (dot(x,x) + df*scale)/rand(Chisq(n+df))
end

#sample variance for other iid random effects
function sampleVCs(mme::MME,sol::Array{Float64,1})
    for effect in  mme.rndTrmVec
        trmi       = effect.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        x          = sol[startPosi:endPosi]
        effect.vcOld  = effect.vcNew
        effect.vcNew  = sample_variance(x,trmi.nLevels, effect.df, effect.scale)
    end
end

#sampel Pi
function samplePi(nEffects, nTotal)
    return rand(Beta(nTotal-nEffects+1, nEffects+1))
end
