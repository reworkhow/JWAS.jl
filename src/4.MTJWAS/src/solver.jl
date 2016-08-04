function Jacobi(A,x,b,p;tolerance=0.000001,output=10)
    D       = diag(A)
    res     = A*x
    resid   = b-res
    tempSol = resid./D
    diff    = sum(resid.^2)
    n    = size(A,1)
    iter = 0
    while ((diff/n > tolerance) & (iter<1000))
        iter += 1
        x = p*tempSol + (1-p)*x
        res     = A*x
        resid   = b-res
        tempSol = resid./D + x
        diff    = sum(resid.^2)
        if iter%output == 0
            println(iter," ",diff/n)
        end
    end
    return x
end

function GaussSeidel(A,x,b;tolerance=0.000001,output=10)
    n = size(x,1)
    for i=1:n
        x[i] = ((b[i] - A[:,i]'x)/A[i,i])[1,1] + x[i]
    end
    diff = sum((A*x-b).^2)
    iter = 0
    while ((diff/n > tolerance) & (iter<1000))
        iter += 1
        for i=1:n
            x[i] = ((b[i] - A[:,i]'x)/A[i,i])[1,1] + x[i]
        end
        diff = sum((A*x-b).^2)
        if iter%output == 0
            println(iter," ",diff/n)
        end
    end
    return x
end

function Gibbs(A,x,b,nIter;outFreq=100) #General Gibbs
    n = size(x,1)
    xMean = zeros(n)
    for iter = 1:nIter
        if iter%outFreq==0
            println("at sample: ",iter)
        end
        for i=1:n
            cVarInv = 1.0/A[i,i]
            cMean   = cVarInv*(b[i] - A[:,i]'x)[1,1] + x[i]
            x[i]    = randn()*sqrt(cVarInv) + cMean
        end
        xMean += (x - xMean)/iter
    end
    return xMean
end

function Gibbs(A,x,b) #one iteration of General Gibbs (NOT \lambda version of MME)
    n = size(x,1)
    for i=1:n
        cVarInv = 1.0/A[i,i]
        cMean   = cVarInv*(b[i] - A[:,i]'x)[1,1] + x[i]
        x[i]    = randn()*sqrt(cVarInv) + cMean
    end
end
