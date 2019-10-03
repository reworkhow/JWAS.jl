function LOOCV(model,adjusted_phenotype,vara,vare)
  X = model.
  X = model.M.genotypes
  y = Array(adjusted_phenotype[:adjusted_phenotype])
  vara,vare=vara,vare

  nind,nmarker=size(model.M.genotypes)
  if nind > nmarker
    res=eMEM(X,y,vara,vare)
  else
    res=eBV(X,y,vara,vare)
  end
  return res
end

#PRESS p<<n
function eMEM(X,y,vara,vare)
    λ = vare/vara
    nInd = size(X,1)
    X = [ones(nInd) X]
    nPar = size(X,2)
    ident = eye(nPar)
    ident[1,1]=0.0
    H0 = inv(X'X+ident*λ)
    β =H0*X'y
    H_array=[(X[i,:]'H0*X[i,:])[1]::Float64 for i=1:nInd]
    e = (y - X*β)./(1-H_array)
    println("mean square of error is ",(e'e/nInd)[1])
    println("prediction accuracy is ",cor(y,y+e))
    return e
end

#PRESS p>>n
function eBV(X,y,vara,vare)
    λ = vare/vara
    nInd,nMarkers = size(X)
    Z = [ones(nInd) eye(nInd)]
    if nInd > nMarkers
        Hsub = eye(nInd)+inv(X*X'+eye(nInd)*0.01)*λ #X*X' is not full rank if n>p
    else
        Hsub = eye(nInd)+inv(X*X')*λ
    end
    H0 = inv([nInd ones(nInd)'
             ones(nInd) Hsub])
    β =H0*Z'y
    H_array=[(Z[i,:]'H0*Z[i,:])[1]::Float64 for i=1:nInd]
    e = (y - Z*β)./(1-H_array)
    println("mean square of error is ",(e'e/nInd)[1])
    println("prediction accuracy is ",cor(y,y+e))
    return e
end

function QInv(X,y,vara,vare)
    bigV = 1e5
    nInd, nMarkers = size(X)
    X = [ones(nInd) X]
    D = diagm([bigV; ones(nMarkers)*vara])
    V = X*D*X'+eye(size(X,1))*vare
    Q =[y'y y'
        y   V];
    Qi = inv(Q)

    SSE = 0.0
    for j=1:length(y)
        ehatj = -Qi[1+j,1]/(Qi[1,1]*Qi[1+j,1+j] - Qi[1,1+j]*Qi[1+j,1])
        SSE+=ehatj^2
    end
    return SSE
end

################################################################################
# Naive LOOCV
################################################################################
#situation n>>p
function MEM(X,y,vara,vare;proportion=1.0) #y_sub=Z_sub*X\beta+e
    λ = vare/vara
    nInd,nMarkers = size(X)
    SSE   = 0.0

    #get a propotion of nInd to make test faster
    #1.0 all; 0.5 half
    nIter = round(Integer,nInd*proportion)
    fixed = ones(nInd-1,1)


    for i in 1:nIter
        sel = fill(true,nInd)
        sel[i] = false

        X_sub = X[sel,:]
        lhs   = [ sum(fixed) fixed'X_sub
                   X_sub'fixed X_sub'X_sub+eye(nMarkers)*λ]
        y_sub = y[sel]
        rhs   = [sum(y_sub)
                 X_sub'y_sub]
        β = inv(lhs)*rhs
        e = y[i] - β[1] - X[i,:]*β[2:end]
        SSE += e[1,1]^2
    end

    return SSE
end

#situation p>>n
function BV(X,y,vara,vare;proportion=1.0) #y_sub=Z_sub*u+e
    λ = vare/vara
    nInd,nMarkers = size(X)
    lhs = eye(nInd)
    SSE = 0.0

    if nInd > nMarkers
        Ginv = inv(X*X'+eye(nInd)*0.01)*λ #X*X' is not full rank if n>p
    else
        Ginv = inv(X*X')*λ
    end
    #get a subset of nInd to make test faster
    nIter=round(Integer,nInd*proportion)
    y_sub= copy(y)
    lhs_sub= copy(lhs)

    fixed = ones(nInd,1)

    for i in 1:nIter
        lhs_sub[i,i]= 0.0
        y_sub[i]=0.0
        fixed[i]=0.0
        lhs = [ nInd-1 fixed'
                  fixed  lhs_sub+Ginv]
        rhs = [sum(y_sub)
                 y_sub]
        u   = inv(lhs)*rhs
        lhs_sub[i,i]=1.0
        y_sub[i]=y[i]
        fixed[i]=1.0

        e = y[i]-u[1]-u[i+1]
        SSE += e[1,1]^2
    end

    return SSE
end
