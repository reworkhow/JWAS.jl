#define samples for WHICH location parameters to output
function outputSamplesFor(mme::MME,trmStr::AbstractString)
    #add model number => "1:age"
    res = []
    for (m,model) = enumerate(mme.modelVec)
        strVec  = split(model,['=','+'])
        strpVec = [strip(i) for i in strVec]
        if trmStr in strpVec
            res = [res;string(m)*":"*trmStr]
        end
    end #"age"->"1:age","2:age"

    for trmStr in res
        trm     = mme.modelTermDict[trmStr]
        samples = MCMCSamples(trm,Array(Float64,1,1))
        push!(mme.outputSamplesVec,samples)
    end
end

function init_sample_arrays(mme::MME,niter)
    #varaince components for residual
    mme.samples4R = zeros(mme.nModels^2,niter)

    #variance components for random polygenic effects
    if mme.ped != 0
        mme.samples4G = zeros(length(mme.pedTrmVec)^2,niter)
    end

    #location parameters for fixed and random effects except markers
    for i in  mme.outputSamplesVec #resize
        trmi = i.term
        i.sampleArray = zeros(trmi.nLevels,niter)
    end

    #variance components for iid random effects
    for i in  mme.rndTrmVec #resize to be size of nTraits
        i.sampleArray = zeros(mme.nModels^2,niter)#Bug maybe many diff
    end
end

#output samples for location parameers
function outputSamples(mme::MME,sol,iter::Int64)
    for i in  mme.outputSamplesVec
        trmi = i.term
        startPosi  = trmi.startPos
        endPosi    = startPosi + trmi.nLevels - 1
        i.sampleArray[:,iter] = sol[startPosi:endPosi]
    end
    for effect in  mme.rndTrmVec
        effect.sampleArray[iter] = effect.vcNew
    end
end


"""
    outputMCMCsamples(mme::MME,trmStr::AbstractString...)

Get samples for specific variables.
"""
function outputMCMCsamples(mme::MME,trmStr::AbstractString...)
    for i in trmStr
      outputSamplesFor(mme,i)
    end
end
