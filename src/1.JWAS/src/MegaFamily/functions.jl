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

function latent_factor_sample_F(mme,df)
    invDf = inv(mme.R)
    # calculate Phi_RE and update F given crrent iteration and Y
    nObs = size(df,1)
    nObsTrait = length(mme.yobs)
    vareRjs = Array{Float64,1}(undef, nObsTrait)
    for obstraitj = 1:length(mme.yobs)
        lambda_j = mme.Λ[:,obstraitj]
        sumlambdaf = sum([lambda_j[k] .* mme.ySparse[(k-1)*nObs+1:k*nObs] for k in 1:length(lambda_j)])
        vareRj = var(mme.yobs[obstraitj] - sumlambdaf)
        vareRjs[obstraitj] = vareRj
    end
    invLhs = invDf + mme.Λ ./vareRjs' * mme.Λ'
    LambdaDYinv = mme.Λ ./vareRjs'
    if mme.M != 0
        for Mi in mme.M
            X2B2F = [Mi.genotypes*Mi.α[traiti] for traiti in 1:mme.K]
            for i in 1:nObs
                mu = invLhs * (LambdaDYinv * getindex.(mme.yobs, nObs) + invDf * getindex.(X2B2F, nObs))
                for k in 1:mme.K
                    mme.ySparse[(k-1)*nObs+i] = mu[k]
                end
            end
        end
    end
    ycorr = mme.ySparse - vec(Matrix(mme.ySparse)-mme.X*mme.sol)
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
