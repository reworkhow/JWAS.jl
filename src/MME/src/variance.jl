function sampleVariance(x, n, df, scale) #for markers
    return (dot(x,x) + df*scale)/rand(Chisq(n+df))
end

function sampleVCs(mme::MME,sol::Array{Float64,1},iter::Int64)
    for effect in  mme.rndTrmVec
        trmi       = effect.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        x          = sol[startPosi:endPosi]
        effect.vcOld  = effect.vcNew
        effect.vcNew  = sampleVariance(x,trmi.nLevels, effect.df, effect.scale)
        effect.sampleArray[iter] = effect.vcNew
    end   
end