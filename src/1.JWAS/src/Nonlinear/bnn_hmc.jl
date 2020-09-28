# in total L hidden layers
# l_l nodes in l-th hidden layer (e.g., l2 nodes in 2nd hidden layer)

###########
# Z0: marker covariate matrix, n by p, each col is a maker

# Z1: latent layer(1), n by l1 matrix, each col is a latent trait
# Z2: latent layer(2), n by l2 matrix, each col is a latent trait
# ...
# ZL: latent layer(L), n by lL matrix, each col is a latent trait

# y: observed trait, vector of n by 1

# Z_all = [Z1,Z2,Z3,...,ZL]  array of matrix
###########

###########
# W0: marker effects, matrix of p by l1, each col is for a latent trait

# W1: weights from latent(1) to latent(2), matrix of l1 by l2, each col is for a latent trait
# W2: weights from latent(2) to latent(3), matrix of l2 by l3, each col is for a latent trait
# ...
# WL: weights from latent(L) to observed trait, vector of lL by 1

# W_all = [W1,...,WL] array of matrix
###########

###########
# Mu1: vector of length l1
# ...
# MuL: vector of length lL
# mu: overall mean of obsered trait

# Mu_all = [Mu1,...,ML] array of array
###########

###########
# Sigma2z1: in hidden layer 1, residual variance of latent trait, vector of length l1
# ...
# Sigma2zL: in hidden layer L, residual variance of latent trait, vector of length lL

# Sigma2z_all = [Sigma2z1,...,Sigma2zL]

# sigma2e: residual variance of observed trait
###########

# nNodes = [3,4,4,3] # nNodes in each hidden layer
# L = length(nNodes) # total number of hidden layer



#helper 1: gradiant of j-th latent trait in l-th hidden layer, zl_j
function calc_gradient_zl_j(j,l,L,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
    # data to be used
    if l==1
        Zl_minus_one=Z0
    else
        Zl_minus_one=Z_all[l-1]
    end

    Zl=Z_all[l]

    if l != L
        Zl_plus_one=Z_all[l+1]
    end

    # weights to be used
    if l==1
        Wl_minus_one = W0
    else
        Wl_minus_one = W_all[l-1]
    end
    Wl = W_all[l]

    # variance to be used
    Sigma2zl          = Sigma2z_all[l]
    if l != L
        Sigma2zl_plus_one = Sigma2z_all[l+1]
    end


    # others to be used
    n  = length(y)

    zl_j=Zl[:,j]
    mul_j = Mu_all[l][j]
    wl_minus_one_j = Wl_minus_one[:,j] #vector of length l_l_minus_one
    if l != L
        Mul_plus_one = mean(Zl_plus_one,dims=1)  #(1,l_l+1)
    end

    sigma2zl_j = Sigma2zl[j]


    # calc_gradient_zl_j = dlogfl_plus_one + dlogfl
    if l==1
        dlogfl = - (zl_j .- mul_j - Zl_minus_one*wl_minus_one_j)/sigma2zl_j
    else
        dlogfl = - (zl_j .- mul_j - tanh.(Zl_minus_one)*wl_minus_one_j)/sigma2zl_j
    end

    if l==L
        wL_j= Wl[j]
        # mu = mean(y)

        dlogfl_plus_one = (y .- mu - tanh.(Zl)*Wl) * (wL_j/sigma2e) .* (-tanh.(zl_j).^2 .+ 1)
    else
        dlogfl_plus_one = (Zl_plus_one - ones(n)*Mul_plus_one - tanh.(Zl)*Wl) * (Wl[j,:] ./ Sigma2zl_plus_one) .* (-tanh.(zl_j).^2 .+ 1)
    end

    gradient_zl_j = dlogfl + dlogfl_plus_one

    return gradient_zl_j  #(n,1)
end


# helper 2: log p(z|y) to help calculate the acceptance rate in sampling of zl_j
function calc_log_p_z(j,l,L,nNodes,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
    # data to be used
    if l==1
        Zl_minus_one=Z0
    else
        Zl_minus_one=Z_all[l-1]
    end

    Zl=Z_all[l]

    if l != L
        Zl_plus_one=Z_all[l+1]
    end


    # weights to be used
    if l==1
        Wl_minus_one = W0
    else
        Wl_minus_one = W_all[l-1]
    end

    Wl = W_all[l]


    # others to be used
    Sigma2zl          = Sigma2z_all[l]
    if l != L
        Sigma2zl_plus_one = Sigma2z_all[l+1]
    end

    n  = length(y)

    zl_j=Zl[:,j]
    mul_j = Mu_all[l][j]
    wl_minus_one_j = Wl_minus_one[:,j] #vector of length l_l_minus_one
    if l != L
        Mul_plus_one = mean(Zl_plus_one,dims=1)  #(1,l_l+1)
    end

    sigma2zl_j = Sigma2zl[j]


    # calc_log_p_z = logfl_plus_one + logfl
    if l==1
        logfl = - 0.5*(zl_j .- mul_j - Zl_minus_one*wl_minus_one_j).^2/sigma2zl_j .- 0.5*log(sigma2zl_j)
    else
        logfl = - 0.5*(zl_j .- mul_j - tanh.(Zl_minus_one)*wl_minus_one_j).^2/sigma2zl_j .- 0.5*log(sigma2zl_j)
    end

    if l==L
        # mu = mean(y)
        logfl_plus_one = -0.5*(y .- mu - tanh.(Zl)*Wl).^2 /sigma2e .- 0.5*log(sigma2e)
    else
        logfl_plus_one = -0.5*(Zl_plus_one - ones(n)*Mul_plus_one - tanh.(Zl)*Wl).^2 * (ones(nNodes[l+1]) ./ Sigma2zl_plus_one) .- 0.5*sum(log.(Sigma2zl_plus_one))
    end

    log_p_z = logfl + logfl_plus_one

    return log_p_z  #(n,1)
end


#helper 3: one iterations of HMC to sample zl_j
function hmc_one_iteration(nLeapfrog,epsilon,j,l,L,nNodes,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
    n = length(y)
    old_Z_all = Z_all
    # z=Z_all[l][:,j]

    #Step1. update phi ~ N(0,M)
    phi = randn(n)  #rand(n,Normal(0,sigma2phi))

    log_p_old = calc_log_p_z(j,l,L,nNodes,Z0,old_Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu) - 0.5*phi.^2  #-0.5*log(sigma2phi)
    #Step2. update (zl_j,phi) from L leapfrog
    phi += 0.5 * epsilon * calc_gradient_zl_j(j,l,L,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
    for leap_i in 1:nLeapfrog
       #(b) full step of theta
       Z_all[l][:,j] += epsilon * phi  # * 1/sigma2phi
       #(c) half step of phi
       if leap_i==nLeapfrog
           phi += 0.5 * epsilon * calc_gradient_zl_j(j,l,L,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
       else
           phi += epsilon * calc_gradient_zl_j(j,l,L,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
       end
    end

    #Step3. acceptance rate
    # Z_all[l][:,j] = z #Z_all with z_new for all indvidual, just to calculate r
    log_p_new = calc_log_p_z(j,l,L,nNodes,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu) - 0.5*phi.^2
    r         = exp.(log_p_new - log_p_old)  # (n,1)

    # r[r.>1] .= 1   # p_jump, for each r, minimum.([r,1]）
    nojump = rand(n) .> r
    Z_all[l][nojump,j] = old_Z_all[l][nojump,j] #Z_all with z_new for jumping indviduals

    return Z_all
end


#helper 4: run multiple HMC inside one Gibbs iteration, to sample zl_j
function hmc_run(j,l,L,nNodes,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu;epsilon=0.1,nLeapfrog=10,nIter_hmc=10)
    for i in 1:nIter_hmc
        Z_all = hmc_one_iteration(nLeapfrog,epsilon,j,l,L,nNodes,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
    end
    return Z_all
end


#helper 5: MME solver to help sample weights & mu
function MME(y,x;lambda=0.00001)
    lhs= x'x + I*lambda
    rhs= x'y
    Ch = cholesky(lhs)
    iL = inv(Ch.L)
    iLhs = inv(Ch)
    p=size(lhs,2)
    sol=iLhs * rhs + iL'randn(p)

    return sol
end

#helper 6: sample weight connected to zl_j in the left (wl_mimus_one_j)
function sample_weights_mu(l,j,L,Z0,Z_all,y,W0,W_all)
    if l==L+1
        y_mme = y
    else
        y_mme = Z_all[l][:,j]
    end

    n     = length(y_mme)

    if l==1
        #sample marker effect w0_j
        x_mme = [ones(n) Z0]
    else
        #sample weights wl_mimus_one_j
        x_mme = [ones(n) tanh.(Z_all[l-1])]
    end

    # lambda = 0.00001 #var(y_mme)/var(wl_minus_one_j) #vary/varw
    mu_weight = MME(y_mme,x_mme)

    mul_j = mu_weight[1]
    wl_minus_one_j = mu_weight[2:end]  #first is mu

    return mul_j,wl_minus_one_j
end


# helper 7: given trained model, and new inputs, calculate predictions
function cal_prediction(Z0,L,W0,W_all,Mu_all,mu)  #Z0 is marker in testing data
    n= size(Z0,1)
    Z= ones(n)*(Mu_all[1]') + Z0*W0  # Z1
    for l in 2:L
        Z = ones(n)*(Mu_all[l]') + tanh.(Z)*W_all[l-1] #Z2...ZL
    end
    y = ones(n)*mu + tanh.(Z)*W_all[L]

    return y
end

# helper 8: given Z1_hat, W_all ,Mu_all, mu -> y_hat
function cal_prediction_fromZ1(Z1,L,W_all,Mu_all,mu)
    n= size(Z1,1)
    Z=Z1

    for l in 2:L
        Z = ones(n)*(Mu_all[l]') + tanh.(Z)*W_all[l-1] #Z2...ZL
    end
    y = ones(n)*mu + tanh.(Z)*W_all[L]

    return y
end

function sample_latent_traits_hmc(yobs,mme,ycorr)  #ycorr is residual
    ylats_old = mme.ySparse         # current values of each latent trait in 1st hidden layer; vec(matrix of n by l1)
    μ_ylats   = mme.ySparse - ycorr # mean of each latent trait in 1st hidden layer; =Z1-residual;  vec(matrix of n by l1)
                                    # = vcat(getEBV(mme,1).+mme.sol[1],getEBV(mme,2).+mme.sol[2]))
    σ2_yobs   = mme.σ2_yobs         # residual variance of yobs (scalar)

    L=mme.L
    nNodes=mme.nNodes

    #reshape the vector to n by l1
    nobs, ntraits = length(mme.obsID), mme.nModels
    ylats_old     = reshape(ylats_old,nobs,ntraits)


    ## already have W0, Mu_all[1]
    ## need to update Z1,W1, Z2,W2, ..., ZL,WL
    ############# START ##################
    # sample Z1
    for j=1:nNodes[1]
        mme.Z_all = hmc_run(j,1,L,nNodes,mme.M[1].genotypes,mme.Z_all,yobs,mme.W0,mme.W_all,mme.Sigma2z_all,σ2_yobs,mme.Mu_all,mme.mu)  #Z0=mme.M[1].genotypes
    end
    # sample W1, Z2,W2, ... ZL-1,WL-1, ZL
    for l=2:L
        #sample latent trait zl_j
        for j=1:nNodes[l]
            ##sample weight and latent trait in l-th hidden layer
            (mme.Mu_all[l][j], mme.W_all[l-1][:,j]) = sample_weights_mu(l,j,L,mme.M[1].genotypes,mme.Z_all,yobs,mme.W0,mme.W_all)
            #sample latent trait zl_j
            mme.Z_all = hmc_run(j,l,mme.L,mme.nNodes,mme.M[1].genotypes,mme.Z_all,yobs,mme.W0,mme.W_all,mme.Sigma2z_all,σ2_yobs,mme.Mu_all,mme.mu)
        end
    end
    #sample WL (l=L)
    (mme.mu, mme.W_all[L]) = sample_weights_mu(L+1,999,L,mme.M[1].genotypes,mme.Z_all,yobs,mme.W0,mme.W_all)

    mme.weights_NN = mme.W_all[L]

    ############# END ##################

    ylats_new        = mme.Z_all[1]

    mme.ySparse = vec(ylats_new)
    ycorr[:]    = mme.ySparse - μ_ylats  # =(ySparse_new - ySparse_old) + ycorr   change of ycorr due to updating Z1

    #sample σ2_yobs (vare)
    yhat = cal_prediction(mme.M[1].genotypes,L,mme.W0,mme.W_all,mme.Mu_all,mme.mu)  #Z0=mme.M[1].genotypes
    residuals = yobs-yhat#[ones(nobs) tanh.(ylats_new)]*weights
    mme.σ2_yobs= dot(residuals,residuals)/rand(Chisq(nobs)) #(dot(x,x) + df*scale)/rand(Chisq(n+df))
end
