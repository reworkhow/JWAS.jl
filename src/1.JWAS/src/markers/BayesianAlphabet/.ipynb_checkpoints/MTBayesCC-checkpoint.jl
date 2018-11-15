
function sampleMarkerEffectsBayesCC!(xArray,xpx,wArray,alphaArray,meanAlphaArray,
                                    deltaArray,meanDeltaArray,
                                    uArray,meanuArray,
                                    invR0,invG0,iIter,BigPi,labels,burnin)

    nMarkers = length(xArray)
    nTraits  = length(alphaArray)
    Ginv     = invG0
    Rinv     = invR0

       α = zeros(nTraits)
    newu = zeros(nTraits)
    oldu = zeros(nTraits)
       δ = zeros(nTraits)
       w = zeros(nTraits) #for rhs


      nlable    = length(labels)
      probDelta = Array(Float64,nlable)
      logDelta = Array(Float64,nlable)
      αlpha     = Array(Array{Float64,1},nlable)
      RinvLhs   = Array(Array{Float64,2},nlable)
      RinvRhs   = Array(Array{Float64,2},nlable)

      for i in 1:length(labels)
        δi = labels[i]
        D  = diagm(δi)
        RinvLhs[i] = D*Rinv*D #split better
        RinvRhs[i] = Rinv*D
      end

    for marker=1:nMarkers

        x    = xArray[marker]

        for trait = 1:nTraits
            α[trait]  = alphaArray[trait][marker]
         oldu[trait]  = newu[trait] = uArray[trait][marker]
            #δ[trait]  = deltaArray[trait][marker]
            w[trait]  = dot(x,wArray[trait])+xpx[marker]*oldu[trait]
        end

        stdnorm = randn(nTraits)
        for label in 1:length(labels)
            lhs       = RinvLhs[label]*xpx[marker]+Ginv
            rhs       = w'*RinvRhs[label]

            #lhsC      = cholfact(lhs)            #if lhs pd
            #lhsC      = chol(lhs)            #if lhs pd
            invLhs    = inv(lhs)                #nTrait X nTrait
            invLhsC   = chol(Hermitian(invLhs))
            #gHat     = lhsC\rhs' #nTrait X 1
            gHat      = invLhs*rhs'
            #probDelta[label]= sqrt(1.0/det(lhsC))*exp(0.5*(rhs*gHat)[1,1])+BigPi[label]
            #logDelta[label]=-0.5*(log(det(lhsC))-(rhs*gHat)[1,1])+log(BigPi[label])
            logDelta[label]=-0.5*(log(det(lhs))-(rhs*gHat)[1,1])+log(BigPi[label])

            αlpha[label]    = gHat + invLhsC'*stdnorm
            #println(lhs,rhs,gHat)
        end

        #println("Before: ",probDelta)

        for label in 1:length(labels)
          deno =0.0
          for i in 1:length(labels)
            deno += exp(logDelta[i]-logDelta[label])
          end
          probDelta[label]=1/deno
        end

        #println("After: ",probDelta)


        #probDelta  = probDelta/sum(probDelta)
        #whichlabel = rand(Categorical(probDelta),Distributions.NoArgCheck())
        whichlabel = rand(Categorical(probDelta))


        δ           = labels[whichlabel]
        α           = αlpha[whichlabel]
        newu        = diagm(δ)*α #α.*δ

        # adjust for locus j
        for trait = 1:nTraits
            BLAS.axpy!(oldu[trait]-newu[trait],x,wArray[trait])
            alphaArray[trait][marker]      = α[trait]
            deltaArray[trait][marker]      = δ[trait]
            uArray[trait][marker]          = newu[trait]
            if iIter>burnin
                meanAlphaArray[trait][marker] += (α[trait] - meanAlphaArray[trait][marker])/(iIter-burnin)
                meanDeltaArray[trait][marker] += (δ[trait] - meanDeltaArray[trait][marker])/(iIter-burnin)
                meanuArray[trait][marker]     += (newu[trait] - meanuArray[trait][marker])/(iIter-burnin)
            end
        end
    end
end
