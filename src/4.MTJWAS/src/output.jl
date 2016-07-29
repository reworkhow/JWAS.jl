#function outputSamplesFor(mme::MME,trmStr::AbstractString)
#    trm  = mme.modelTermDict[trmStr]
#    samples = MCMCSamples(trm,Array(Float64,1,1))
#    push!(mme.outputSamplesVec,samples)
#end

function init_sample_arrays(mme::MME,niter)
    samples4R = Array(Array{Float64,2},niter)
    samples4G  = []
    if mme.ped != 0
        samples4G = Array(Array{Float64,2},niter)
    end
    MCMCsamples(samples4R,samples4G)
    #for i in  mme.outputSamplesVec
    #    trmi = i.term
    #    i.sampleArray = zeros(niter,trmi.nLevels)
    #end
    #for i in  mme.rndTrmVec
    #    trmi = i.term
    #    i.sampleArray = zeros(niter)
    #end
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
