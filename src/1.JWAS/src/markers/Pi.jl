################################################################################
#
# Sample parameter for the indicator variable : π
#
################################################################################
#single-trait
function samplePi(nEffects::Number, nTotal::Number)
    return rand(Beta(nTotal-nEffects+1, nEffects+1))
end

#multi-trait
function samplePi(deltaArray::Union{Array{Array{Float64,1},1},Array{Array{Float32,1},1}},BigPi)
  temp = deltaArray[1]
  ntraits = size(deltaArray,1)
  for traiti = 2:ntraits
    temp = [temp deltaArray[traiti]]
  end

  iloci = 1
  nLoci_array=zeros(length(BigPi))
  for i in keys(BigPi) #assume order of key won't change
    temp2 = broadcast(-,temp,i')
    nLoci =  sum(mean(abs.(temp2),dims=2).==0.0)
    nLoci_array[iloci] = nLoci +1
    iloci = iloci +1
  end
  tempPi = rand(Dirichlet(nLoci_array))

  iloci = 1
  for i in keys(BigPi)
    BigPi[i] = tempPi[iloci]
    iloci = iloci +1
  end
end
