"""
    solve(mme::MME,df::DataFrame;solver="default",printout_frequency=100,tolerance = 0.000001,maxiter = 5000)

* Solve the mixed model equations (no marker information) without estimating variance components.
Available solvers include `default`, `Jacobi`, `Gauss-Seidel`, and `Gibbs sampler`.
"""
function solve(mme::MME,
               df::DataFrame;
               solver="default",
               printout_frequency=100,
               tolerance = 0.000001,
               maxiter = 5000,
               heterogeneous_residuals = false,
               double_precision=false)
    df = check_pedigree_genotypes_phenotypes(mme,df,heterogeneous_residuals)
    set_default_priors_for_variance_components(mme,df)
    if size(mme.mmeRhs)==()
        getMME(mme,df)
    end
    if double_precision == true
        A = map(Float64,mme.mmeLhs)
        b = map(Float64,mme.mmeRhs)
        x = zeros(Float64,size(mme.mmeRhs,1))
    else
        A = map(Float32,mme.mmeLhs)
        b = map(Float32,mme.mmeRhs)
        x = zeros(Float32,size(mme.mmeRhs,1))
    end
    if solver=="Jacobi"
        return [getNames(mme) Jacobi(A,x,b,
                                    tolerance=tolerance,
                                    maxiter=maxiter,
                                    printout_frequency=printout_frequency)]
    elseif solver=="Gauss-Seidel"
        return [getNames(mme) GaussSeidel(A,x,b,
                              tolerance=tolerance,
                              maxiter=maxiter,
                              printout_frequency=printout_frequency)]
    elseif solver=="Gibbs" && mme.nModels !=1
        return [getNames(mme) Gibbs(A,x,b,
                              maxiter,printout_frequency=printout_frequency)]
    elseif solver=="Gibbs" && mme.nModels==1
        return [getNames(mme) Gibbs(A,x,b,mme.R,
                              maxiter,printout_frequency=printout_frequency)]
    elseif solver=="default"
        println("To solve the equations, please choose a solver. (run ?solve for help)")
        println("Following values are returned:")
        println("(1) names ;(2) incidence matrix;")
        println("(3) left-hand side and (4) right-hand side of mixed model equations.")
        return [getNames(mme),mme.X,A,b]
    else
        error("Please try solver=`default`,`Jacobi`,`Gauss-Seidel`, or `Gibbs sampler`\n")
    end
end

################################################################################
#Solvers including Jacobi, Gauss-Seidel, and Gibbs (general,lambda,one iteration)
################################################################################
function Jacobi(A,x,b,p=0.7;tolerance=0.000001,printout_frequency=10,maxiter=1000)
    n       = size(A,1)   #number of linear equations
    D       = diag(A)
    error   = b - A*x
    diff    = sum(error.^2)/n

    iter    = 0
    while (diff > tolerance) && (iter < maxiter)
        iter   += 1
        error   = b - A*x
        x_temp  = error./D + x
        x       = p*x_temp + (1-p)*x
        diff    = sum(error.^2)/n

        if iter%printout_frequency == 0
            println("at iteration ",iter,": ",diff)
        end
    end
    return x
end

function GaussSeidel(A,x,b;tolerance=0.000001,printout_frequency=10,maxiter=1000)
    n = size(A,1)
    for i = 1:n
        x[i] = (b[i] - A[:,i]'x)/A[i,i] + x[i]
    end
    error = b - A*x
    diff  = sum(error.^2)/n

    iter  = 0
    while (diff > tolerance) & (iter < maxiter)
        iter += 1
        for i = 1:n
            x[i] = (b[i] - A[:,i]'x)/A[i,i] + x[i]
        end

        error = b - A*x
        diff  = sum(error.^2)/n
        if iter%printout_frequency == 0
            println("at iteration ",iter,": ",diff)
        end
    end
    return x
end

#Gibbs for lambda version of MME (single-trait)
function Gibbs(A,x,b,vare::AbstractFloat,niter::Int64;printout_frequency=100)
    n = size(x,1)
    xmean = zeros(n)
    for iter = 1:niter
        if iter%printout_frequency==0
            println("at iteration: ",iter)
        end
        for i=1:n
            invlhs = 1.0/A[i,i]
            μ      = invlhs*(b[i] - A[:,i]'x) + x[i]
            x[i]   = randn()*sqrt(invlhs*vare) + μ
        end
        xmean += (x - xmean)/iter
    end
    return xmean
end

#General Gibbs (multi-trait)
function Gibbs(A,x,b,niter::Int64;printout_frequency=100)
    n = size(x,1)
    xmean = zeros(n)
    for iter = 1:niter
        if iter%printout_frequency==0
            println("at iteration: ",iter)
        end
        for i=1:n
            if A[i,i] != 0.0 #issue70, zero diagonals in MME
                invlhs = 1.0/A[i,i]
                μ      = invlhs*(b[i] - A[:,i]'x) + x[i]
                x[i]   = randn()*sqrt(invlhs) + μ
            end
        end
        xmean += (x - xmean)/iter
    end
    return xmean
end

#one iteration of Gibbs sampler for lambda version of MME (single-trait)
function Gibbs(A,x,b,vare::AbstractFloat)
    for i = 1:size(x,1)
        if A[i,i] != 0.0 #issue70, zero diagonals in MME
            invlhs  = 1.0/A[i,i]
            μ       = invlhs*(b[i] - A[:,i]'x) + x[i]
            x[i]    = randn()*sqrt(invlhs*vare) + μ
        end
    end
end

#one iteration of Gibbs sampler for general version of MME (multi-trait)
function Gibbs(A,x,b)
    for i = 1:size(x,1)
        if A[i,i] != 0.0 #issue70, zero diagonals in MME
            invlhs  = 1.0/A[i,i]
            μ       = invlhs*(b[i] - A[:,i]'x) + x[i]
            x[i]    = randn()*sqrt(invlhs) + μ
        end
    end
end
