function summary(X::Array{Array{Float64,2},1};index=false)
  Xmean = mean(X)
  Xvar  = mean(X.^2)-Xmean.^2

  if index != false
    xArray = zeros(length(X))
    i=1
    for thisX in X
      xArray[i]=thisX[index[1],index[2]]
      i+=1
    end
    println(describe(xArray))
    plt[:hist](xArray);
  end
  return Dict("mean"=>Xmean,"variance"=>Xvar)
end

function summary(X::Array{Array{Float64,1},1})
  Xmean = mean(X)
  Xvar  = mean(X.^2)-Xmean.^2
  return Dict("mean"=>Xmean,"variance"=>Xvar)
end


#export summary
