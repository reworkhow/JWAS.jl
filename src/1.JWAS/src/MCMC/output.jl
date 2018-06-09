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


"""
    get_outputX_others(model)

Make incidence matrices for effects involve in EBV inclung J, ϵ, pedTrmVec except marker covariates
"""
function get_outputX_others(model)
    #trick to avoid errors (PedModule.getIDs(ped) [nongeno ID;geno ID])
    model.output_X["ϵ"]=mkmat_incidence_factor(model.output_ID,
                                               PedModule.getIDs(model.ped))[:,1:length(model.ped.setNG)]
    #model.output_X["J"] is convenient in SSBRrun
    for i in model.pedTrmVec
        model.output_X[i]=mkmat_incidence_factor(model.output_ID,
                                                model.modelTermDict[i].names)
    end
end

# function getEBV(model)
#     for traiti in 1:mme.nModels
#
#         res["Posterior mean of marker effects"]
#
#        file= output_file*"_"*"marker_effects_"*string(mme.lhsVec[traiti])*".txt"
#        res["EBV"*"_"*string(mme.lhsVec[traiti])]=
#           misc.get_breeding_values(mme,file,header=(mme.M.markerID[1]!="NA")?true:false)[1]
#     end
#     output_others(mme)
# end
