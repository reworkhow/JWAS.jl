function outputSamplesFor(mme::MME,trmStr::AbstractString)
    trm  = mme.modelTermDict[trmStr]
    samples = MCMCSamples(trm,Array(Float64,1,1))
    push!(mme.outputSamplesVec,samples)
end

function initSampleArrays(mme::MME,niter)
    mme.resVarSampleArray = zeros(niter)
    if mme.ped != 0
        n = size(mme.GiNew,1)^2
        mme.genVarSampleArray = zeros(niter,n)
    end
    for i in  mme.outputSamplesVec
        trmi = i.term
        i.sampleArray = zeros(niter,trmi.nLevels)
    end
    for i in  mme.rndTrmVec
        trmi = i.term
        i.sampleArray = zeros(niter)
    end  
end

#output samples for location parameers
function outputSamples(mme::MME,sol,iter::Int64)
    for i in  mme.outputSamplesVec
        trmi = i.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        i.sampleArray[iter,:] = sol[startPosi:endPosi]
    end
end