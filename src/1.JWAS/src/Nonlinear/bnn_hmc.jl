#tianjing: modified from non-linear.jl
#sample all hidden nodes together

#                 |---- Z[:,1] ----- Z0*W0[:,1]
#yobs ---f(X)----|---- Z[:,2] ----- Z0*W0[:,2]
#                 |---- Z[:,3] ----- Z0*W0[:,3]
#
# in total 1 hidden layers, with l1 nodes.

###########
# X: marker covariate matrix, n by p, each col is for 1 marker
# Z: latent traits, n by l1 matrix, each col is for 1 latent trait
# y: observed trait, vector of n by 1
###########

###########
# W0: marker effects, matrix of p by l1, each col is for a latent trait
# W1: weights from hidden layer to observed trait, vector of l1 by 1
###########

###########
# Mu0: bias of latent traits, vector of length l1
# mu: bias of observed trait, scaler
###########

###########
# Sigma2z: residual variance of latent trait, diagonal matrix of size l1*l1
# sigma2e: residual variance of observed trait, scaler
###########


#helper 1: calculate gradiant of all latent traits for all individual
function calc_gradient_z(Z,y,W0,W1,Sigma2z,sigma2e,Mu0,mu,ycorr)
    n=length(y)
    tanhZ=tanh.(Z)
    #dlogfz =- (Z - ones(n)*Mu0' - X*W0) * inv(Sigma2z)  #(n,l1)
    dlogfz = - ycorr * inv(Sigma2z)
    dlogfy = ((y .- mu - tanhZ*W1)/sigma2e) * W1' .* (-tanhZ.^2 .+ 1) #(n,l1)
    gradient_z = dlogfz + dlogfy

    return gradient_z  #(n,l1)
end


# helper 2: calculate log p(z|y) to help calculate the acceptance rate
function calc_log_p_z(Z,y,W0,W1,Sigma2z,sigma2e,Mu0,mu,ycorr)

    n=length(y)
    #logfz = -0.5*sum(((Z-ones(n)*Mu0'-X*W0).^2)*inv(Sigma2z),dims=2) .- (0.5*log(prod(diag(Sigma2z))))
    logfz = -0.5*sum((ycorr.^2)*inv(Sigma2z),dims=2) .- (0.5*log(prod(diag(Sigma2z))))  #
    logfy = -0.5*(y .- mu - tanh.(Z)*W1).^2 /sigma2e .- 0.5*log(sigma2e)

    log_p_z = logfz + logfy

    return log_p_z  #(n,1)
end


#helper 3: one iterations of HMC to sample Z
function hmc_one_iteration(nLeapfrog,epsilon,Z,y,W0,mu_W1,Sigma2z,sigma2e,Mu0,ycorr)
    mu=mu_W1[1]
    W1=mu_W1[2:end]
    n,l1 = size(Z)
    old_Z = Matrix(Z)

    #Step1. update phi ~ N(0,M)
    phi = randn(n,l1)  #rand(n,Normal(0,sigma2phi))
    log_p_old = calc_log_p_z(old_Z,y,W0,W1,Sigma2z,sigma2e,Mu0,mu,ycorr) - 0.5*sum(phi.^2,dims=2)  #(n,1)
    #Step2. update (zl,phi) from 10 leapfrog
    phi += 0.5 * epsilon * calc_gradient_z(Z,y,W0,W1,Sigma2z,sigma2e,Mu0,mu,ycorr)  #(n,l1)
    for leap_i in 1:nLeapfrog
       #(b) full step of theta
       Z += epsilon * phi  # (n,l1)
       ycorr += epsilon * phi #update ycorr due to change of Z
       #(c) half step of phi
       if leap_i==nLeapfrog
           phi += 0.5 * epsilon * calc_gradient_z(Z,y,W0,W1,Sigma2z,sigma2e,Mu0,mu,ycorr)
       else
           phi += epsilon * calc_gradient_z(Z,y,W0,W1,Sigma2z,sigma2e,Mu0,mu,ycorr)
       end
    end

    #Step3. acceptance rate
    log_p_new = calc_log_p_z(Z,y,W0,W1,Sigma2z,sigma2e,Mu0,mu,ycorr) - 0.5*sum(phi.^2,dims=2) #(n,1)
    r         = exp.(log_p_new - log_p_old)  # (n,1)
    nojump = rand(n) .> r  # bool (n,1)

    for i in 1:n
        if nojump[i]
            Z[i,:] = old_Z[i,:]
        end
    end

    return Z
end



function sample_latent_traits_hmc(yobs,mme,ycorr,iter)  #ycorr is residual
    # ylats_old = mme.ySparse         # current values of each latent trait; vec(matrix of n by l1)
    # μ_ylats   = mme.ySparse - ycorr # mean of each latent trait; =Z-residual;  vec(matrix of n by l1)
    #                                  # = vcat(getEBV(mme,1).+mme.sol[1],getEBV(mme,2).+mme.sol[2]))
    # σ2_yobs   = mme.σ2_yobs         # residual varianum_latent_traits=mme.M[1].ntraits

    #reshape the vector to n by l1
    nobs, ntraits = length(mme.obsID), mme.nModels
    ylats_old     = reshape(ylats_old,nobs,ntraits)
    ycorr         = reshape(ycorr,nobs,ntraits)

    ############# START ##################
    # sample latent trait (Z)
    ylats_new = hmc_one_iteration(10,0.1,mme.Z,yobs,mme.W0,mme.weights_NN,mme.R,σ2_yobs,mme.Mu0,ycorr)

    #sample weights (W1)
    M = [ones(nobs) tanh.(mme.Z)]
    lhs= M'M + I*0.00001   #here +I*0.00001 is just to keep lhs numerical stable
    rhs= M'yobs
    Ch = cholesky(lhs)
    iL = inv(Ch.L)
    iLhs = inv(Ch)
    weights=iLhs * rhs + iL'randn(size(M,2))*sqrt(σ2_yobs)  #construct a vector alpha_hat+(L')^{-1}z, z~N(0,vare)
    mme.weights_NN = weights

    ############# END ##################
    mme.ySparse = vec(ylats_new)
    ycorr[:]    = mme.ySparse - μ_ylats  # =(ySparse_new - ySparse_old) + ycorr   change of ycorr due to updating Z

    #sample σ2_yobs (vare)
    residuals = yobs .- mme.mu - tanh.(ylats_new)*mme.W1
    σ2_yobs= dot(residuals,residuals)/rand(Chisq(nobs)) #(dot(x,x) + df*scale)/rand(Chisq(n+df))
    mme.σ2_yobs= σ2_yobs
    mme.vare_mean += (σ2_yobs - mme.vare_mean)/iter

end
