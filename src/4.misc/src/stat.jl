"""
    report(X::Array{Array{Float64,2},1};index=false)
* show summary statistics for MCMC samples (matrices or index [i,j] of matrices)
"""
function report(X::Array{Array{Float64,2},1};index=false)
  if index == false
    Xmean = mean(X)
    Xvar  = mean(X.^2)-Xmean.^2
    println("Summary Stats:")
    println("Mean:\n",round(Xmean,digits=6))
    println("Variance:\n",round(Xvar,digits=6))
  else
    xArray = zeros(length(X))
    i=1
    for thisX in X
      xArray[i]=thisX[index[1],index[2]]
      i+=1
    end
    println(DataFrames.describe(xArray))
    #plt[:hist](xArray,normed=true,bins=round(Int,length(xArray)/20));
    #Plots.histogram(xArray,nbins=round(Int,length(xArray)/20),normalize=true)
    #Plots.title!("Trait $(index[1]) and Trait $(index[2])")
  end
end

"""
    report(X::Array{Array{Float64,1},1};index=false)
* show summary statistics for MCMC samples (vectors or index i of vectors)
"""
function report(X::Array{Array{Float64,1},1};index=false)
  if index==false
    println("Summary Stats:")
    println("Mean:\n",round(mean(X),digits=6))
  else
    xArray = zeros(length(X))
    i=1
    for thisX in X
      xArray[i]=thisX[index]
      i+=1
    end
    println(DataFrames.describe(xArray))
    #plt[:hist](xArray,normed=true,bins=round(Int,length(xArray)/20));
    #Plots.histogram(xArray,nbins=round(Int,length(xArray)/20),normalize=true)
    #Plots.title!("Trait $index")
  end
end

export report
