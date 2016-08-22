function ssBayesC0_constantvariance(matrices::HybridMatrices,
                 geno::misc.Genotypes,fixed::misc.FixedMatrix,
                 ped::PedModule.Pedigree,
                 input::misc.InputParameters;outFreq=5000)

    y   = matrices.y.full
    X   = matrices.X.full
    W   = matrices.W.full
    Zn  = matrices.Z.n
    Ai_nn = matrices.Ai.nn

    current   = misc.Current(input,geno,fixed,y)
    current.fixed_effects = zeros(length(current.fixed_effects)+1) #add one for μ_g
    current.imputation_residual = zeros(matrices.num.pedn)

    output                 = misc.Output(input,geno,fixed)
    output.mean_fixed_effects = zeros(length(output.mean_fixed_effects)+1) #add one for μ_g
    output.mean_imputation_residual        = zeros(matrices.num.pedn)

    wGibbs = misc.GibbsMats(W)
    xGibbs = misc.GibbsMats(X)

    #if input.estimateVariance == false
      #construct lhs for sampleEpsilon!
      λ_ϵ             = input.varResidual/input.varGenotypic
      Z_n             = matrices.Z._n
      lhs_ϵ           = Z_n'Z_n+Ai_nn*λ_ϵ
      lhsCol_ϵ        = [lhs_ϵ[:,i] for i=1:size(lhs_ϵ,1)]
      lhsDi_ϵ         = 1.0./diag(lhs_ϵ)
      sd_ϵ            = sqrt(lhsDi_ϵ*current.varResidual)
      #construct lhs for sample marker effects
      λ               = current.varResidual/current.varEffect
      lhsDi           = [1.0/(wGibbs.xpx[i]+λ)::Float64 for i=1:size(wGibbs.xpx,1)]
      sd              = sqrt(lhsDi*current.varResidual)
    #end

    println("running ",input.method," with a MCMC of length ",input.chainLength)

    for iter = 1:input.chainLength

      current.iter += 1
      # sample fixed effects
      misc.sample_fixed!(xGibbs,current,output)
      # sample marker effects
      misc.sample_random_ycorr!(wGibbs,current,output,lhsDi,sd)
      # sample epsilon
      SSBR.sampleEpsi!(matrices,current,output,lhsCol_ϵ,lhsDi_ϵ,sd_ϵ)

      if (iter%outFreq ==0)
          println("This is iteration ",iter)
      end
    end

    estimatedMarkerEffects = output.meanMarkerEffects

    mu_g = output.mean_fixed_effects[end]
    EBV = matrices.J.full*mu_g+matrices.M.full*estimatedMarkerEffects
    EBV[1:matrices.num.pedn,:] += output.mean_imputation_residual

    IDs=PedModule.getIDs(ped);
    EBV=DataFrame(ID=IDs,EBV=vec(EBV))

    return Results(matrices,output,EBV)
end

export ssGibbs
