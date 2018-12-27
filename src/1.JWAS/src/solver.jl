"""
    solve(mme::MME,df::DataFrame;solver="default",printout_frequency=100,tolerance = 0.000001,maxiter = 5000)

* Solve the mixed model equations (no marker information) without estimating variance components.
Available solvers include `default`, `Jacobi`, `GaussSeidel`, and `Gibbs sampler`.
"""
function solve(mme::MME,
               df::DataFrame;
               solver="default",
               printout_frequency=100,
               tolerance = 0.000001,
               maxiter = 5000)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    p = size(mme.mmeRhs,1)
    if solver=="Jacobi"
        return [getNames(mme) Jacobi(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,0.3,
                                    tolerance=tolerance,maxiter=maxiter,
                                    outFreq=printout_frequency)]
    elseif solver=="GaussSeidel"
        return [getNames(mme) GaussSeidel(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,
                              tolerance=tolerance,maxiter=maxiter,
                              outFreq=printout_frequency)]
    elseif solver=="Gibbs" && mme.nModels !=1
        return [getNames(mme) Gibbs(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,
                              maxiter,outFreq=printout_frequency)]
    elseif solver=="Gibbs" && mme.nModels==1
        return [getNames(mme) Gibbs(mme.mmeLhs,fill(0.0,p),mme.mmeRhs,mme.RNew,
                              maxiter,outFreq=printout_frequency)]
    elseif solver=="default"
        return [getNames(mme) mme.mmeLhs\mme.mmeRhs]
    else
        error("No this solver. Please try `default`,`Jacobi`,`GaussSeidel`, or `Gibbs sampler`\n")
    end
end

################################################################################
#Solvers including Jacobi, GaussSeidel, and Gibbs (general,lambda,one_iteration)
################################################################################
function Jacobi(A,x,b,p;tolerance=0.000001,outFreq=10,maxiter=1000)
    D       = diag(A)
    res     = A*x
    resid   = b-res
    tempSol = resid./D
    diff    = sum(resid.^2)
    n       = size(A,1)
    iter    = 0
    while ((diff/n > tolerance) & (iter<maxiter))
        iter   += 1
        x       = p*tempSol + (1-p)*x
        res     = A*x
        resid   = b-res
        tempSol = resid./D + x
        diff    = sum(resid.^2)
        if iter%outFreq == 0
            println(iter," ",diff/n)
        end
    end
    return x
end

function GaussSeidel(A,x,b;tolerance=0.000001,outFreq=10,maxiter=1000)
    n = size(x,1)
    for i=1:n
        x[i] = (b[i] - A[:,i]'x)/A[i,i] + x[i]
    end
    diff = sum((A*x-b).^2)
    iter = 0
    while ((diff/n > tolerance) & (iter<maxiter))
        iter += 1
        for i=1:n
            x[i] = (b[i] - A[:,i]'x)/A[i,i] + x[i]
        end
        diff = sum((A*x-b).^2)
        if iter%outFreq == 0
            println(iter," ",diff/n)
        end
    end
    return x
end

#Gibbs for \lambda version of MME (single-trait)
function Gibbs(A,x,b,varRes::Float64,nIter::Int64;outFreq=100)
    n = size(x,1)
    xMean = zeros(n)
    for iter = 1:nIter
        if iter%outFreq==0
            println("at sample: ",iter)
        end
        for i=1:n
            cVarInv = 1.0/A[i,i]
            cMean   = cVarInv*(b[i] - A[:,i]'x) + x[i]
            x[i]    = randn()*sqrt(cVarInv*varRes) + cMean
        end
        xMean += (x - xMean)/iter
    end
    return xMean
end

#General Gibbs (multi-trait)
function Gibbs(A,x,b,nIter::Int64;outFreq=100)
    n = size(x,1)
    xMean = zeros(n)
    for iter = 1:nIter
        if iter%outFreq==0
            println("at sample: ",iter)
        end
        for i=1:n
            cVarInv = 1.0/A[i,i]
            cMean   = cVarInv*(b[i] - A[:,i]'x) + x[i]
            x[i]    = randn()*sqrt(cVarInv) + cMean
        end
        xMean += (x - xMean)/iter
    end
    return xMean
end

#one iteration of Gibbs for \lambda version of MME (single-trait)
function Gibbs(A,x,b,varRes::Float64)
    n = size(x,1)
     for i=1:n
        cVarInv = 1.0/A[i,i]
        cMean   = cVarInv*(b[i] - A[:,i]'x) + x[i]
        x[i]    = randn()*sqrt(cVarInv*varRes) + cMean
    end
end

#one iteration of Gibbs for general version of MME (multi-trait)
function Gibbs(A,x,b)
    n = size(x,1)
    for i=1:n
      if A[i,i] != 0 # get rid of it #double-check
        cVarInv = 1.0/A[i,i]
        cMean   = cVarInv*(b[i] - A[:,i]'x) + x[i]
        x[i]    = randn()*sqrt(cVarInv) + cMean
      end
    end
end
