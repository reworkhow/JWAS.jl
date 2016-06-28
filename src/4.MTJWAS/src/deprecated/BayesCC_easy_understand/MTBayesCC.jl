
function sampleMarkerEffectsBayesCC!(xArray,xpx,wArray,alphaArray,meanAlphaArray,
                                    deltaArray,meanDeltaArray,
                                    uArray,meanuArray,
                                    invR0,invG0,iIter,BigPi)

    nMarkers = length(xArray)
    nTraits  = length(alphaArray)
    Ginv     = invG0
    Rinv     = invR0

       α = zeros(nTraits)
    newu = zeros(nTraits)
    oldu = zeros(nTraits)
       δ = zeros(nTraits)

       w = zeros(nTraits) #for rhs

      Delta  = Dict{Array{Float64,1},Float64}()
      LhsRhs = Dict{Array{Float64,1},Array{Array{Float64,2},1}}()

    for j=1:nMarkers

        x    = xArray[j]

        for k=1:nTraits
               α[k] = alphaArray[k][j]
            oldu[k] = newu[k] = uArray[k][j]
        #       δ[k] = deltaArray[k][j]
        end
        #oldu = copy(newu) #move in

        for trait = 1:nTraits
            w[trait] = dot(x,wArray[trait])+xpx[j]*oldu[trait]
        end

        for δ in keys(BigPi)
            #lhs       = δ.*Rinv.*δ'*xpx[j]+Ginv
            #rhs       = w'*Rinv.*δ'
            k1        = diagm(δ)
            k2        = Rinv*k1
            lhs       = k1*k2*xpx[j]+Ginv
            rhs       = w'*k2


            lhsC      = cholfact(lhs)            #if lhs pd
            invLhs    = inv(lhsC)                #nTrait X nTrait

            gHat       = invLhs*rhs' #nTrait X 1
            Delta[δ]  = sqrt(1.0/det(lhsC))*exp(0.5*(rhs*gHat)[1,1])+BigPi[δ]
            LhsRhs[δ] = Array[gHat,invLhs] #expensive, copy
        end
        probDelta = cumsum(collect(values(Delta)))
        probDelta = probDelta/probDelta[end]

        #whichδ      = sum(rand() .> probDelta)+1 ##****
        whichδ      = findfirst(x->x>rand(),probDelta)

        δ           = copy(collect(keys(Delta))[whichδ]) ##copy required
                                                         ##otherwise modifiying keys(Delta)

        gHat,invLhs = LhsRhs[δ]
        L           = chol(invLhs)
        α           = gHat + L*rand(nTraits)

        #=
        k1        = diagm(δ)
        k2        = Rinv*k1
        lhs       = k1*k2*xpx[j]+Ginv
        rhs       = w'*k2
        lhsC      = cholfact(lhs)
        gHat      = inv(lhsC)*rhs'
        L         = inv(lhsC[:U])
        α         = gHat + L*randn(nTraits)
        =#

        newu        = diagm(δ)*α #α.*δ

        # adjust for locus j
        for trait = 1:nTraits
            BLAS.axpy!(oldu[trait]-newu[trait],x,wArray[trait])
            meanAlphaArray[trait][j] += (α[trait] - meanAlphaArray[trait][j])/iIter
            meanDeltaArray[trait][j] += (δ[trait] - meanDeltaArray[trait][j])/iIter
            meanuArray[trait][j]     += (newu[trait] - meanuArray[trait][j])/iIter

            alphaArray[trait][j]      = α[trait]
            deltaArray[trait][j]      = δ[trait]
            uArray[trait][j]          = newu[trait]
        end
    end
end
