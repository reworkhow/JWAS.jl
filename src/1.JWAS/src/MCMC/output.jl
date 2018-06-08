################################################################################
#User-interface to output EBV
################################################################################
"""
    outputEBV(model,IDs::Array{String,1})

Output estimated breeding values and prediction error variances for IDs.
"""
function outputEBV(model,IDs::Array{String,1})
    println("Estimated breeding values and prediction error variances will be included in the output.")
    model.output_ID=IDs
end


function output_others(model)
    ϵnames = model.modelTermDict["1:ϵ"].names
    #trick to avoid errors (model.output_ID NOT ⊂ ϵ)
    model.output_X["ϵ"]=mkmat_incidence_factor(model.output_ID,[ϵnames;model.output_ID])[:,1:length(ϵnames)]
    #model.output_X["J"] is convenient in SSBRrun
    for i in model.pedTrmVec
        model.output_X[i]=mkmat_incidence_factor(model.output_ID,
                                                model.modelTermDict[i].names)
    end
end
