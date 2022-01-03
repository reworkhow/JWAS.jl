#Below function is to re-phase model for latent factors (based on function from Tianjing)
function factor_model_equation(model_equations,K)
    # model equations: y = intercept + geno

    lhs, rhs = strip.(split(model_equations,"="))
    model_equations = ""
    # new equations: f1=intercept+geno; f2=intercept+geno

    for i = 1:K
        model_equations = model_equations*"f"*string(i)*"="*rhs*";"
    end
    model_equations = model_equations[1:(end-1)] #remove the semicolon at the end
end

function latent_factors_setup!(mme)
    Lambda = mme.Lamb
    t = length(mme.Lamb.trait_names)
    K = length(Lambda.factor_names)
    if Lambda.λ == false
        Lambda.λ = [zeros(t) for i in 1:K]
        Lambda.is_estimate = true
    end
    Lambda.β = [copy(Lambda.λ[factori]) for factori = 1:K]
    Lambda.γ = [ones(typeof(Lambda.λ[factori][1]),t) for factori = 1:K]
    Lambda.π = zeros(K)
    Lambda.δ = [1;zeros(K-1)]

    Lambda.ycorr_obsTrait = copy(mme.yobs)
    Lambda.meanLambda  = [zero(Lambda.λ[factori]) for factori = 1:K]
    Lambda.meanLambda2 = [zero(Lambda.λ[factori]) for factori = 1:K]
    Lambda.meanGamma   = [zero(Lambda.γ[factori]) for factori = 1:K]
    Lambda.mean_pi            = zeros(K)
    Lambda.mean_pi2           = zeros(K)
    Lambda.meanDelta          = zeros(K)
    Lambda.meanDelta2         = zeros(K)
end

function latent_factor_sample_F(mme,df,ycorr)
    invDf = inv(mme.R)
    # calculate Phi_RE and update F given crrent iteration and Y
    nObs = length(mme.obsID)
    Lambda = mme.Lamb
    nObsTrait = length(Lambda.trait_names)

    Λ = vcat([x' for x in Lambda.λ]...)
    ycorr[:] = ycorr + mme.X*mme.sol # fk
    for obstraitj = 1:nObsTrait
        lambda_j = Λ[:,obstraitj]
        sumlambdaf = sum([lambda_j[k] .* ycorr[(k-1)*nObs+1:k*nObs] for k in 1:length(lambda_j)])
        Lambda.ycorr_obsTrait[obstraitj] = mme.yobs[obstraitj] - sumlambdaf # used to calculate varRjs
    end
    invLhs = invDf + Λ ./Lambda.varRjs' * Λ'
    LambdaDYinv = Λ ./Lambda.varRjs'
    if mme.M != 0
        for Mi in mme.M
            X2B2F = [Mi.genotypes*Mi.α[traiti] for traiti in 1:Lambda.K]
            for i in 1:nObs
                mu = invLhs * (LambdaDYinv * getindex.(mme.yobs, i) + invDf * getindex.(X2B2F, i))
                samples = rand(MvNormal(mu, convert(Array,Symmetric(invLhs))))
                for k in 1:Lambda.K
                    ycorr[(k-1)*nObs+i] = samples[k]
                end
            end
        end
    end
    ycorr[:] =  ycorr-mme.X*mme.sol

    if mme.M != 0
        for Mi in mme.M
            for traiti in 1:Mi.ntraits
                if Mi.α[traiti] != zero(Mi.α[traiti])
                    ycorr[(traiti-1)*Mi.nObs+1 : traiti*Mi.nObs] =
                    ycorr[(traiti-1)*Mi.nObs+1 : traiti*Mi.nObs] - Mi.genotypes*Mi.α[traiti]
                end
            end
        end
    end
    return ycorr
end


function BayesC_Lambda_sampler!(fArray,yCorr,λ,β,γ,τ,varRjs,π)
    logPi         = log(π)
    logPiComp     = log(1-π)
    logGamma0     = logPi
    invVarRjs     = 1 ./varRjs
    varEffects    = 1/τ .* varRjs
    invVarEffects = 1 ./  varEffects
    logVarEffects = log.(varEffects)
    nElements     = length(λ)

    fpf = dot(fArray,fArray)
    for j=1:nElements
        rhs = (dot(fArray,yCorr[j]) + fpf*λ[j])*invVarRjs[j]
        lhs = fpf*invVarRjs[j] + invVarEffects[j]
        invLhs = 1/lhs
        gHat   = rhs*invLhs
        logGamma1  = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp ## Why?
        probGamma1 = 1/(1+ exp(logGamma0 - logGamma1))
        oldLambda  = λ[j]

        if(rand()<probGamma1)
            γ[j] = 1
            β[j] = gHat + randn()*sqrt(invLhs)
            λ[j] = β[j]
            BLAS.axpy!(oldLambda-λ[j],fArray,yCorr[j])
        else
            if (oldLambda!=0)
                BLAS.axpy!(oldLambda,fArray,yCorr[j])
            end
            γ[j] = 0
            β[j] = randn()*sqrt(varEffects[j])
            λ[j] = 0
        end
    end
end

function BayesC_Lambda_sampler_parallel!(ycorr,nobs,Lambda)
    Threads.@threads for k in 1:Lambda.K #nfactors
         BayesC_Lambda_sampler!(ycorr[(k-1)*nobs+1:k*nobs],Lambda.ycorr_obsTrait,
                                Lambda.λ[k], Lambda.β[k], Lambda.γ[k], prod(Lambda.δ[1:k]),
                                Lambda.varRjs, Lambda.π[k])
    end
end


function delta_sampler!(Lambda)
    t = length(Lambda.trait_names)
    for l in 2:Lambda.K
        delta_a = Lambda.delta_a + t*(Lambda.K-l+1)/2
        delta_b = Lambda.delta_b
        for k in l:Lambda.K
            delta_h = Lambda.δ[1:k]
            delta_h_corr = prod(deleteat!(delta_h, l))
            sumlambda2 = sum(Lambda.λ[k].^2 ./ Lambda.varRjs)
            delta_b += 0.5*delta_h_corr * sumlambda2
        end
        Lambda.δ[l] = rand(Gamma(delta_a, 1/delta_b))
    end
end

function rescaleF!(ycorr, Lambda, mme)
    nObs = length(mme.obsID)
    ycorr2 = ycorr.^2
    F_sizes = [mean(ycorr2[(k-1)*nObs+1:k*nObs]) for k in 1:Lambda.K]
    for k in 1:Lambda.K
    ycorr[(k-1)*nObs+1:k*nObs] = ycorr[(k-1)*nObs+1:k*nObs] ./ sqrt(F_sizes[k])
    end
    Mi = mme.M[1]
    for k in 1:Lambda.K
        Mi.α[k] = Mi.α[k] ./ sqrt(F_sizes[k])
    end
    mme.R = mme.R ./ F_sizes'
    delta_factor = [F_sizes[1];exp.(diff(log.(F_sizes)))]
    Lambda.δ = delta_factor
end
