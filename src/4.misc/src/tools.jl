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

function center2(X) #without 2 here, not use center! and center togehter, maybe Julia bug
    X=copy(X)
    nrow,ncol = size(X)
    colMeans = mean(X,1)
    BLAS.axpy!(-1,ones(nrow)*colMeans,X)
    return colMeans
end

#borrow from julia class

function mkmat_incidence_factor(b)
    factor=unique(b)
    coMat= spzeros(length(b),length(factor))

    dictFactor = Dict()
    index=1
    for i in factor
        dictFactor[i]=index
        index+=1
    end

    nrow=1
    for i in b
        myindex=dictFactor[i]
        coMat[nrow,myindex]=1
        nrow=nrow+1
    end
    return full(coMat),factor
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

function mkmat_incidence_factor(b)
    factor=unique(b)
    coMat= spzeros(length(b),length(factor))

    dictFactor = Dict()
    index=1
    for i in factor
        dictFactor[i]=index
        index+=1
    end

    nrow=1
    for i in b
        myindex=dictFactor[i]
        coMat[nrow,myindex]=1
        nrow=nrow+1
    end
    return full(coMat),factor
end

function delete_rc(x,rowi,colj)
    numRow,numCol = size(x)
    nrow = [1:numRow]
    ncol = [1:numCol]

    if rowi != 0
        deleteat!(nrow,rowi)
    end

    if colj != 0
        deleteat!(ncol,colj)
    end

    return x[nrow,ncol]
end

function find_column_byname(file,names)
  colname = split(readline(file))
  nameIndex=Int64[]

  for i in names
    index=find(colname,names)
    push!(nameIndex,index)
  end

  return nameIndex
end

function get_gene_freq(X,id=true)
  ncol=size(X,2)
  start=1
  freq=[]

  if id==true
    freq=zeros(size(X,2)-1)
    start += 1
  else
    freq=zeros(size(X,2))
  end

  k=1
  for i=start:ncol
    freq[k]=mean(X[:,i])/2
    k+=1
  end

  return freq
end

function getXpRinvX(X, Rinv)
    ncol = size(X)[2]
    XpRinvX = [((X[:,i].*Rinv)'X[:,i])[1]::Float64 for i=1:ncol]
    return XpRinvX
end

function getXpRinvX(X)
    XpRinvX = [dot(X[:,i],X[:,i]) for i=1:size(X,2)]
    return XpRinvX
end


