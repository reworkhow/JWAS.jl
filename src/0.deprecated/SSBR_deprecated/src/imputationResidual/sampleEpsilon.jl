function sampleEpsi!(mats::HybridMatrices,current::misc.Current,out::misc.Output)
    yCorr         = current.yCorr
    varRes        = current.varResidual
    λ             = current.varResidual/current.varGenotypic
    iIter         = 1/current.iter
    Z_n           = mats.Z._n
    Ai_nn         = mats.Ai.nn
    meanEpsi      = out.mean_imputation_residual
    ϵ             = current.imputation_residual

    #add back {Zn'Zn}_{i,i} *ϵ
    current.yCorr = current.yCorr + Z_n*ϵ
    rhs = Z_n'current.yCorr

    #println("varGenotypic is ",current.varGenotypic)
    #println("varRes is ",current.varResidual)

    lhs = Z_n'Z_n+Ai_nn*λ #better to construct Z_n'Z_n outside

    sample_random_rhs!(lhs,rhs,current,out) #use this general function for sample epsilon(Gianola Book)

    current.yCorr = current.yCorr - Z_n*ϵ
end

#with known variances
function sampleEpsi!(mats::HybridMatrices,current::misc.Current,out::misc.Output,lhsCol,lhsDi,sd)
    Z_n  = mats.Z._n
    yCorr= current.yCorr
    ϵ    = current.imputation_residual


    #maybe current.yCorr = current.yCorr + Z_n*ϵ better
    yCorr[:] = yCorr[:] + Z_n*ϵ
    #rhs = Z_n'current.yCorr
    rhs = Z_n'yCorr

    sample_random_rhs!(lhsCol,rhs,current,out,lhsDi,sd)

    #current.yCorr = current.yCorr - Z_n*ϵ
    yCorr[:] = yCorr[:] - Z_n*ϵ
end
