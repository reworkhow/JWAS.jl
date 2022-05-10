###############################################################################
#MT-BayesCC is required for some scenarios, e.g, only [0,0] and [1,1] for delta
#in this case, MT-BayesC won't work.
#***** this one is still needed for research projects
###############################################################################
function sampleMarkerEffectsBayesCC!(xArray,xpx,wArray,alphaArray,
                                    deltaArray,
                                    uArray,
                                    invR0,invG0,iIter,BigPi,labels)

    nMarkers = length(xArray)
    ntraits  = length(alphaArray)
    Ginv     = invG0
    Rinv     = invR0

       α = zeros(ntraits)
    newu = zeros(ntraits)
    oldu = zeros(ntraits)
       δ = zeros(ntraits)
       w = zeros(ntraits) #for rhs


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

        for trait = 1:ntraits
            α[trait]  = alphaArray[trait][marker]
         oldu[trait]  = newu[trait] = uArray[trait][marker]
            #δ[trait]  = deltaArray[trait][marker]
            w[trait]  = dot(x,wArray[trait])+xpx[marker]*oldu[trait]
        end

        stdnorm = randn(ntraits)
        for label in 1:length(labels)
            lhs       = RinvLhs[label]*xpx[marker]+Ginv
            rhs       = w'*RinvRhs[label]
            invLhs    = inv(lhs)                #ntraits X ntraits
            invLhsC   = chol(Hermitian(invLhs))
            #gHat     = lhsC\rhs' #ntraits X 1
            gHat      = invLhs*rhs'
            #probDelta[label]= sqrt(1.0/det(lhsC))*exp(0.5*(rhs*gHat)[1,1])+BigPi[label]
            #logDelta[label]=-0.5*(log(det(lhsC))-(rhs*gHat)[1,1])+log(BigPi[label])
            logDelta[label]=-0.5*(log(det(lhs))-(rhs*gHat)[1,1])+log(BigPi[label])

            αlpha[label]    = gHat + invLhsC'*stdnorm
        end


        for label in 1:length(labels)
          deno =0.0
          for i in 1:length(labels)
            deno += exp(logDelta[i]-logDelta[label])
          end
          probDelta[label]=1/deno
        end

        #probDelta  = probDelta/sum(probDelta)
        #whichlabel = rand(Categorical(probDelta),Distributions.NoArgCheck())
        whichlabel = rand(Categorical(probDelta))


        δ           = labels[whichlabel]
        α           = αlpha[whichlabel]
        newu        = diagm(δ)*α #α.*δ

        # adjust for locus j
        for trait = 1:ntraits
            BLAS.axpy!(oldu[trait]-newu[trait],x,wArray[trait])
            alphaArray[trait][marker]      = α[trait]
            deltaArray[trait][marker]      = δ[trait]
            uArray[trait][marker]          = newu[trait]
        end
    end
end

function samplePi(deltaArray,BigPi,BigPiMean,iter,labels)
  temp = deltaArray[1]
  ntraits = size(deltaArray,1)
  for traiti = 2:ntraits
    temp = [temp deltaArray[traiti]]
  end

  iloci = 1
  nLoci_array=zeros(BigPi)
  for i in labels
    temp2 = broadcast(-,temp,i')
    nLoci =  sum(mean(abs(temp2),2).==0.0)
    nLoci_array[iloci] = nLoci +1
    iloci = iloci +1
  end

  BigPi[:] = rand(Dirichlet(nLoci_array))
  BigPiMean[:] += (BigPi-BigPiMean)/iter
end

function setPi(Pi)
  #label with the index of Array
  #labels[1]=[0.0,0.0],labels[2]=[1.0,1.0]
  #BigPi[1]=0.3,BigPi[2]=0.7
  nlabels= length(Pi)
  labels = Array(Array{Float64,1},nlabels)
  BigPi  = Array(Float64,nlabels)
  whichlabel=1
  for pair in sort(collect(Pi), by=x->x[2],rev=true)
    key=pair[1]
    labels[whichlabel]=copy(key)
    BigPi[whichlabel]=copy(pair[2])
    whichlabel = whichlabel+1
  end
  BigPiMean = zeros(BigPi)
  return labels,BigPi,BigPiMean
end
