################################################################################
#User-interface to set training and testing individuals
################################################################################
#function set_training(model,IDs)
#    IDs = map(String,vec(IDs)) #Array{String,1}
#    model.training_ID=IDs
#end


"""
    outputEBV(model,IDs::Array;PEV=false)

Output estimated breeding values and prediction error variances (defaulting to false) for IDs.
"""
function outputEBV(model,IDs;PEV=false)
    #print("Estimated breeding values and prediction error variances will be included in the output.")
    IDs = map(String,vec(IDs)) #Array{String,1}
    model.output_ID=IDs
end

function set_testing(model,IDs;PEV=false)
    #print("Estimated breeding values and prediction error variances will be included in the output.")
    IDs = map(String,vec(IDs)) #Array{String,1}
    model.output_ID=IDs
end

"""
    get_outputX_others(model)

Make incidence matrices for effects involve in EBV inclung J, ϵ, pedTrmVec except marker covariates
"""
function get_outputX_others(model,single_step_analysis)
    #trick to avoid errors (PedModule.getIDs(ped) [nongeno ID;geno ID])
    if single_step_analysis == true
        model.output_X["ϵ"]=mkmat_incidence_factor(model.output_ID,
                            PedModule.getIDs(model.ped))[:,1:length(model.ped.setNG)]
        #model.output_X["J"] is in SSBRrun
    end
    if model.pedTrmVec != 0
        for i in model.pedTrmVec
            model.output_X[i]=mkmat_incidence_factor(model.output_ID,
                                                     model.modelTermDict[i].names)
        end
    end
end

export outputEBV
