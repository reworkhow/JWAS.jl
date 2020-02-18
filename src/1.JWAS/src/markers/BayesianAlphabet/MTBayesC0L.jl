function MTBayesL!(xArray,xRinvArray,xpRinvx,
                   wArray,alphaArray,gammaArray,
                   vare,varEffect)
    invR0    = inv(vare)
    invG0    = inv(varEffect)
    nMarkers = length(xArray)
    nTraits  = length(alphaArray)
    Lhs      = zeros(nTraits,nTraits)
    Rhs      = zeros(nTraits)
    newα     = zeros(typeof(alphaArray[1][1]),nTraits)
    oldα     = zeros(typeof(alphaArray[1][1]),nTraits)


    function getGi_function(gammaArray)
        f1(x,j,Gi)  = Gi ./ x[j]
        f2(x,j,Gi)  = Gi
        size(gammaArray,1)>1 ? f1 : f2
    end
    getGi = getGi_function(gammaArray)

    for marker=1:nMarkers
        x, xRinv = xArray[marker], xRinvArray[marker]
        for trait=1:nTraits #unadjust for locus j
            oldα[trait] = newα[trait] = alphaArray[trait][marker]
            Rhs[trait]  = dot(xRinv,wArray[trait])+xpRinvx[marker]*oldα[trait]
        end
        Rhs = invR0*Rhs
        Lhs = xpRinvx[marker]*invR0 + getGi(gammaArray,marker,invG0)
        for trait = 1:nTraits
            lhs  = Lhs[trait,trait]
            ilhs = 1/lhs
            rhs  = Rhs[trait] - (Lhs[trait,:]'newα)[1]
            mu   = ilhs*rhs + newα[trait]
            alphaArray[trait][marker] = mu + randn()*sqrt(ilhs)
            newα[trait] = alphaArray[trait][marker]
        end
        for trait = 1:nTraits
            BLAS.axpy!(oldα[trait]-newα[trait],x,wArray[trait])
        end
    end
end

MTBayesC0!(xArray,xRinvArray,xpRinvx,wArray,alphaArray,vare,varEffect) =
  MTBayesL!(xArray,xRinvArray,xpRinvx,wArray,alphaArray,[1.0],vare,varEffect)

function sampleGammaArray!(gammaArray,alphaArray,mmeMG)
    Gi = inv(mmeMG)
    nMarkers = size(gammaArray,1)
    nTraits  = length(alphaArray[1])==1 ? 1 : length(alphaArray)

    Q  = zeros(nMarkers)
    nTraits > 1 ? calcMTQ!(Q,nMarkers,nTraits,alphaArray,Gi) : calcSTQ!(Q,nMarkers,alphaArray,Gi)
    gammaDist = Gamma(0.5,4) # 4 is the scale parameter, which corresponds to a rate parameter of 1/4
    candidateArray = 1 ./ rand(gammaDist,nMarkers)
    uniformArray = rand(nMarkers)
    acceptProbArray = exp.(Q ./4 .*(2 ./ gammaArray - candidateArray))
    replace = uniformArray .< acceptProbArray
    gammaArray[replace] = 2 ./ candidateArray[replace]
end

function calcMTQ!(Q,nMarkers,nTraits,alphaArray,Gi)
    for locus = 1:nMarkers
      for traiti = 1:nTraits
          for traitj = 1:nTraits
              Q[locus] += alphaArray[traiti][locus]*alphaArray[traitj][locus]*Gi[traiti,traitj]
          end
      end
    end
end

function calcSTQ!(Q,nMarkers,alphaArray,Gi)
    Q .= alphaArray.^2 ./Gi
end
