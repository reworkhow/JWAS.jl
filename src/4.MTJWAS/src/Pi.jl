#Convert Genetic variance to marker effect variance based on Pi
function genetic2marker(M::Genotypes,Pi::Dict)
  G=M.G #genetic variance
  nTraits = size(G,1)
  denom   = zeros(nTraits,nTraits)
  for i in 1:nTraits
    for j in i:nTraits
      pi_selected = filter((k,v)->k[i]==1.0 && k[j]==1.0,Pi)
      denom[i,j] = M.sum2pq*sum(values(pi_selected))
      denom[j,i] = denom[i,j]
    end
  end
  M.G = M.G ./ denom
  M.G_is_marker_variance = true
end


function samplePi(deltaArray,BigPi,BigPiMean,iter)
  temp = deltaArray[1]
  nTraits = size(deltaArray,1)
  for traiti = 2:nTraits
    temp = [temp deltaArray[traiti]]
  end

  iloci = 1
  nLoci_array=zeros(2^nTraits)
  for i in keys(BigPi) #assume order of key won't change
    temp2 = broadcast(-,temp,i')
    nLoci =  sum(mean(abs(temp2),2).==0.0)
    nLoci_array[iloci] = nLoci +1
    iloci = iloci +1
  end
  tempPi = rand(Dirichlet(nLoci_array))

  iloci = 1
  for i in keys(BigPi)
    BigPi[i] = tempPi[iloci]
    BigPiMean[i] += (tempPi[iloci]-BigPiMean[i])/iter
    iloci = iloci +1
  end
end

function samplePi(deltaArray,BigPi,BigPiMean,iter,labels)
  temp = deltaArray[1]
  nTraits = size(deltaArray,1)
  for traiti = 2:nTraits
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

function samplePi(nEffects, nTotal)#for single trait
    return rand(Beta(nTotal-nEffects+1, nEffects+1))
end
