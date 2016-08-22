function samplePi(nEffects, nTotal)
    return rand(Beta(nTotal-nEffects+1, nEffects+1))
end

function sample_variance(x, n, df, scale)
    return (dot(x,x) + df*scale)/rand(Chisq(n+df))
end

function sample_epsilon_variance(ϵ,Ai,n,df,scale) #assume df same to vareffect
    return ((ϵ'Ai*ϵ)[1,1] + df*scale)/rand(Chisq(n+df))
end

function sample_fixed!(mats::GibbsMats,current::Current,out::Output)
    α             = current.fixed_effects
    meanAlpha     = out.mean_fixed_effects

    nObs,nEffects = mats.nrows,mats.ncols
    xArray        = mats.xArray
    xpx           = mats.xpx
    yCorr         = current.yCorr
    varRes        = current.varResidual
    iIter         = 1/current.iter

    for j=1:nEffects
        x = xArray[j]
        rhs = dot(x,yCorr) + xpx[j]*α[j,1]
        lhs      = xpx[j]
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j,1]
        α[j]     = mean + randn()*sqrt(invLhs*varRes)
        BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
        meanAlpha[j] += (α[j] - meanAlpha[j])*iIter
    end
end

function sample_random_rhs!(lhs,rhs,current::Current,out::Output) #(Gianola Book)
    α             = current.imputation_residual
    meanAlpha     = out.mean_imputation_residual

    nEffects      = size(lhs,1)
    varRes        = current.varResidual
    iIter         = 1/current.iter

    for i = 1:nEffects #argument lhs here is a sparse matrix for whole lhs
        α[i] = 0.0
        rhsi = rhs[i] - lhs[:,i]'α #column-major
        lhsi = lhs[i,i]
        invLhs = 1.0/lhsi
        meani  = invLhs*rhsi[1]
        #println(invLhs*varRes)
        #if invLhs*varRes < 0
        #  println(lhsi,"  ",varRes)
        #end
        α[i] = meani + randn()*sqrt(invLhs*varRes)
        meanAlpha[i] += (α[i] - meanAlpha[i])*iIter
    end
end

function sample_random_ycorr!(mats::GibbsMats,current::Current,out::Output)#sample vare and vara
    α             = current.α
    meanAlpha     = out.meanMarkerEffects

    nObs,nEffects = mats.nrows,mats.ncols
    xArray        = mats.xArray
    xpx           = mats.xpx
    yCorr         = current.yCorr
    varRes        = current.varResidual
    λ             = current.varResidual/current.varEffect
    iIter         = 1/current.iter

    for j=1:nEffects
        x = xArray[j]
        rhs = dot(x,yCorr) + xpx[j]*α[j,1]
        lhs      = xpx[j] + λ
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j,1]
        α[j]     = mean + randn()*sqrt(invLhs*varRes)
        BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
        meanAlpha[j] += (α[j] - meanAlpha[j])*iIter
    end
end




## NOT sample vare and vara
## construct lhsDi,sd outside

function sample_random_ycorr!(mats::GibbsMats,current::Current,out::Output,lhsDi,sd)#not sample vare and vara
    α             = current.α
    meanAlpha     = out.meanMarkerEffects #not good, not general
                                          #Put option in arguments
    nObs,nEffects = mats.nrows,mats.ncols
    xArray        = mats.xArray
    xpx           = mats.xpx
    yCorr         = current.yCorr
    iIter         = 1/current.iter

    for j=1:nEffects
        @inbounds x = xArray[j]
        @inbounds rhs = dot(x,yCorr) + xpx[j]*α[j,1]
        @inbounds mean     = lhsDi[j]*rhs
        @inbounds oldAlpha = α[j,1]
        @inbounds α[j]     = mean + randn()*sd[j]
        @inbounds BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
        @inbounds meanAlpha[j] += (α[j] - meanAlpha[j])*iIter
    end
end

function sample_random_rhs!(lhsCol,rhs,current::Current,out::Output,lhsDi,sd)
    α             = current.imputation_residual
    meanAlpha     = out.mean_imputation_residual

    nEffects      = length(lhsCol)      #arguments lhs here is a array of cols of lhs
    iIter         = 1/current.iter

    for i = 1:nEffects
        @inbounds α[i] = 0.0
        @inbounds rhsi = rhs[i] - lhsCol[i]'α
        @inbounds α[i] = lhsDi[i]*rhsi[1] + randn()*sd[i]
        @inbounds meanAlpha[i] += (α[i] - meanAlpha[i])*iIter
    end
end

export sample_fixed!,sample_random_rhs!,sample_random_ycorr!,sample_variance
export sample_epsilon_variance,samplePi
