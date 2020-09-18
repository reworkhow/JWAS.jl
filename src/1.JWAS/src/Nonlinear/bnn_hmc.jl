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
    z=Z_all[l][:,j]

    #Step1. update phi ~ N(0,M)
    phi = randn(n)  #rand(n,Normal(0,sigma2phi))
    log_p_old = calc_log_p_z(j,l,L,nNodes,Z0,old_Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu) - 0.5*phi.^2  #-0.5*log(sigma2phi)
    #Step2. update (zl_j,phi) from L leapfrog
    phi += 0.5 * epsilon * calc_gradient_zl_j(j,l,L,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
    for leap_i in 1:nLeapfrog
       #(b) full step of theta
       z += epsilon * phi  # * 1/sigma2phi
       #(c) half step of phi
       if leap_i==nLeapfrog
           phi += 0.5 * epsilon * calc_gradient_zl_j(j,l,L,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
       else
           phi += epsilon * calc_gradient_zl_j(j,l,L,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
       end
    end

    #Step3. acceptance rate
    Z_all[l][:,j] = z #Z_all with z_new for all indvidual, just to calculate r
    log_p_new = calc_log_p_z(j,l,L,nNodes,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu) - 0.5*phi.^2
    r         = exp.(log_p_new - log_p_old)  # (n,1)

    # r[r.>1] .= 1   # p_jump, for each r, minimum.([r,1]）
    jump = rand(n) .< r
    old_Z_all[l][jump,j] = z[jump] #Z_all with z_new for jumping indviduals

    return old_Z_all
end


#helper 4: run multiple HMC inside one Gibbs iteration, to sample zl_j
function hmc_run(j,l,L,nNodes,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu;epsilon=0.1,nLeapfrog=10,nIter_hmc=10)
    for i in 1:nIter_hmc
        Z_all = hmc_one_iteration(nLeapfrog,epsilon,j,l,L,nNodes,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
    end
    return Z_all
end


# one iteration in Gibbs sampling: sample all weights and latent traits
function Gibbs_one_iteration(L,nNodes,Z0,Z_all,y,W0,W_all,Mu_all,mu,Sigma2z_all,sigma2e)
    for l=1:L
        ##sample weight and latent trait in l-th hidden layer
        for j=1:nNodes[l]
            #sample weight connected to zl_j in the left (wl_mimus_one_j)
            if l==1
                #sample marker effect w0_j
                (Mu_all[1][j], W0[:,j]) = sample_weights_mu(l,j,L,Z0,Z_all,y,W0,W_all)
            else
                #sample weights wl_mimus_one_j
                (Mu_all[l][j], W_all[l-1][:,j]) = sample_weights_mu(l,j,L,Z0,Z_all,y,W0,W_all)
            end

            #sample latent trait zl_j
            Z_all = hmc_run(j,l,L,nNodes,Z0,Z_all,y,W0,W_all,Sigma2z_all,sigma2e,Mu_all,mu)
        end
    end

    #sample WL (l=L)
    (mu, W_all[L]) = sample_weights_mu(L+1,999,L,Z0,Z_all,y,W0,W_all)

    return Z_all,W0,W_all,mu,Mu_all
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
