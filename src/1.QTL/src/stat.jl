function describe(X::Array{Array{Float64,2},1};index=false)
  if index == false
    Xmean = mean(X)
    Xvar  = mean(X.^2)-Xmean.^2
    println("Summary Stats:")
    println("Mean:\n",round(Xmean,6))
    println("Variance:\n",round(Xvar,6))
  else
    xArray = zeros(length(X))
    i=1
    for thisX in X
      xArray[i]=thisX[index[1],index[2]]
      i+=1
    end
    println(DataFrames.describe(xArray))
    plt[:hist](xArray,normed=true,bins=round(Int,length(xArray)/20));
  end
end

function describe(X::Array{Array{Float64,1},1};index=false)
  if index==false
    println("Summary Stats:")
    println("Mean:\n",round(mean(X),6))
  else
    xArray = zeros(length(X))
    i=1
    for thisX in X
      xArray[i]=thisX[index]
      i+=1
    end
    println(DataFrames.describe(xArray))
    plt[:hist](xArray,normed=true,bins=round(Int,length(xArray)/20));
  end
end


#export summary
