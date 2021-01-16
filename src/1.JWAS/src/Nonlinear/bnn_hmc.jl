#                 |---- Z1[:,1] ----- Z0*W0[:,1]
#yobs ---f(Z1)----|---- Z1[:,2] ----- Z0*W0[:,2]
#                 |---- Z1[:,3] ----- Z0*W0[:,3]
#
# in total 1 hidden layers, with l1 nodes.

###########
# Z0: marker covariate matrix, n by p, each col is for 1 marker
# Z1: latent traits, n by l1 matrix, each col is for 1 latent trait
# y: observed trait, vector of n by 1
###########

###########
# W0: marker effects, matrix of p by l1, each col is for a latent trait
# W1: weights from hidden layer to observed trait, vector of l1 by 1
###########

###########
# Mu1: bias of latent trait, vector of length l1
# mu: bias of observed trait, scaler
###########

###########
# Sigma2z1: residual variance of latent trait, vector of length l1
# sigma2e: residual variance of observed trait, scaler
###########


#helper 1: calculate gradiant of j-th latent trait
function calc_gradient_z1_j(j,Z0,Z1,y,W0,W1,Sigma2z1,sigma2e,Mu1,mu)
    z1_j=Z1[:,j]
    mu1_j = Mu1[j]
    w0_j = W0[:,j] #vector of length p
    sigma2z1_j = Sigma2z1[j]

    dlogfz1 = - (z1_j .- mu1_j - Z0*w0_j)/sigma2z1_j
    w1_j= W1[j]

    dlogfy = (y .- mu - tanh.(Z1)*W1) * (w1_j/sigma2e) .* (-tanh.(z1_j).^2 .+ 1)

    gradient_z1_j = dlogfz1 + dlogfy

    return gradient_z1_j  #(n,1)
end


# helper 2: calculate log p(z|y) to help calculate the acceptance rate
function calc_log_p_z(j,Z0,Z1,y,W0,W1,Sigma2z1,sigma2e,Mu1,mu)
    z1_j=Z1[:,j]
    mu1_j = Mu1[j]
    w0_j = W0[:,j] #vector of length p
    sigma2z1_j = Sigma2z1[j]

    logfz1 = - 0.5*(z1_j .- mu1_j - Z0*w0_j).^2/sigma2z1_j .- 0.5*log(sigma2z1_j)
    logfy = -0.5*(y .- mu - tanh.(Z1)*W1).^2 /sigma2e .- 0.5*log(sigma2e)

    log_p_z1 = logfz1 + logfy

    return log_p_z1  #(n,1)
end


#helper 3: one iterations of HMC to sample z1_j
function hmc_one_iteration(nLeapfrog,epsilon,j,Z0,Z1,y,W0,W1,Sigma2z1,sigma2e,Mu1,mu)
    n = length(y)
    old_Z1 = Z1

    #Step1. update phi ~ N(0,M)
    phi = randn(n)  #rand(n,Normal(0,sigma2phi))
    log_p_old = calc_log_p_z(j,Z0,old_Z1,y,W0,W1,Sigma2z1,sigma2e,Mu1,mu) - 0.5*phi.^2  #-0.5*log(sigma2phi)
    #Step2. update (zl_j,phi) from 10 leapfrog
    phi += 0.5 * epsilon * calc_gradient_z1_j(j,Z0,Z1,y,W0,W1,Sigma2z1,sigma2e,Mu1,mu)
    for leap_i in 1:nLeapfrog
       #(b) full step of theta
       Z1[:,j] += epsilon * phi  # * 1/sigma2phi
       #(c) half step of phi
       if leap_i==nLeapfrog
           phi += 0.5 * epsilon * calc_gradient_z1_j(j,Z0,Z1,y,W0,W1,Sigma2z1,sigma2e,Mu1,mu)
       else
           phi += epsilon * calc_gradient_z1_j(j,Z0,Z1,y,W0,W1,Sigma2z1,sigma2e,Mu1,mu)
       end
    end

    #Step3. acceptance rate
    log_p_new = calc_log_p_z(j,Z0,Z1,y,W0,W1,Sigma2z1,sigma2e,Mu1,mu) - 0.5*phi.^2
    r         = exp.(log_p_new - log_p_old)  # (n,1)

    nojump = rand(n) .> r
    Z1[nojump,j] = old_Z1[nojump,j] #Z1with z_new for jumping indviduals

    return Z1
end



function sample_latent_traits_hmc(yobs,mme,ycorr,iter)  #ycorr is residual
    ylats_old = mme.ySparse         # current values of each latent trait; vec(matrix of n by l1)
    μ_ylats   = mme.ySparse - ycorr # mean of each latent trait; =Z1-residual;  vec(matrix of n by l1)
                                    # = vcat(getEBV(mme,1).+mme.sol[1],getEBV(mme,2).+mme.sol[2]))
    σ2_yobs   = mme.σ2_yobs         # residual variance of yobs (scalar)
    num_latent_traits=mme.M[1].ntraits
    Z0=mme.M[1].genotypes
    Sigma2z1=diag(mme.R)

    #reshape the vector to n by l1
    nobs, ntraits = length(mme.obsID), mme.nModels
    ylats_old     = reshape(ylats_old,nobs,ntraits)


    ############# START ##################
    # sample latent trait (Z1)
    for j=1:num_latent_traits
        mme.Z1 = hmc_one_iteration(10,0.1,j,Z0,mme.Z1,yobs,mme.W0,mme.W1,Sigma2z1,σ2_yobs,mme.Mu1,mme.mu)
    end

    #sample weights (W1)
    X = [ones(nobs) tanh.(mme.Z1)]
    lhs= X'X + I*0.00001   #here +I*0.00001 is just to keep lhs numerical stable
    rhs= X'yobs
    Ch = cholesky(lhs)
    iL = inv(Ch.L)
    iLhs = inv(Ch)
    mu_weight=iLhs * rhs + iL'randn(size(X,2))*sqrt(σ2_yobs)  #construct a vector alpha_hat+(L')^{-1}z, z~N(0,vare)

    mme.mu = mu_weight[1]
    mme.W1 = mu_weight[2:end]
    mme.weights_NN = mme.W1

    ############# END ##################

    ylats_new   = mme.Z1
    mme.ySparse = vec(ylats_new)
    ycorr[:]    = mme.ySparse - μ_ylats  # =(ySparse_new - ySparse_old) + ycorr   change of ycorr due to updating Z1

    #sample σ2_yobs (vare)
    residuals = yobs .- mme.mu - tanh.(ylats_new)*mme.W1
    σ2_yobs= dot(residuals,residuals)/rand(Chisq(nobs)) #(dot(x,x) + df*scale)/rand(Chisq(n+df))
    mme.σ2_yobs= σ2_yobs
    mme.vare_mean += (σ2_yobs - mme.vare_mean)/iter

end
