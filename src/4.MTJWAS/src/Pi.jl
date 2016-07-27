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
