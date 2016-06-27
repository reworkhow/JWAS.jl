type MarkerMatrix
    X::Array{Float64,2}
    xArray::Array{Array{Float64,1},1}
    markerMeans::Array{Float64,2}
    mean2pq::Float64
    G::Float64
    centered::Bool
    function MarkerMatrix(X::Array{Float64,2},G::Float64)
        markerMeans = center!(X) #centering
        p           = markerMeans/2.0
        mean2pq     = (2*p*(1-p)')[1,1]
        xArray      = get_column_ref(X)
        new(X,xArray,markerMeans,mean2pq,G,true)
    end
end

function addMarkers(mme::MME,df,G::Float64)
    M = convert(Array,df)
    mme.M = MarkerMatrix(M,G)
end

function sample_effects_ycorr!(X,xArray,xpx,yCorr,α,meanAlpha,vRes,vEff,iIter)#sample vare and vara
    nObs,nEffects = size(X)
    λ    = vRes/vEff
    for j=1:nEffects
        x        = xArray[j]
        rhs      = dot(x,yCorr) + xpx[j]*α[j,1]
        lhs      = xpx[j] + λ
        invLhs   = 1.0/lhs
        mean     = invLhs*rhs
        oldAlpha = α[j,1]
        α[j]     = mean + randn()*sqrt(invLhs*vRes)
        BLAS.axpy!(oldAlpha-α[j,1],x,yCorr)
        meanAlpha[j] += (α[j] - meanAlpha[j])/iIter
    end
end

function get_column(X,j)
    nrow,ncol = size(X)
    if j>ncol||j<0
        error("column number is wrong!")
    end
    indx = 1 + (j-1)*nrow
    ptr = pointer(X,indx)
    pointer_to_array(ptr,nrow)
end

function get_column_ref(X)
    ncol = size(X)[2]
    xArray = Array(Array{Float64,1},ncol)
    for i=1:ncol
        xArray[i] = get_column(X,i)
    end
    return xArray
end

function center!(X)
    nrow,ncol = size(X)
    colMeans = mean(X,1)
    BLAS.axpy!(-1,ones(nrow)*colMeans,X)
    return colMeans
end