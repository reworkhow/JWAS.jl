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
function calc_gradient_z(ylats,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)
    μ1, w1     = weights_NN[1], weights_NN[2:end]
    g_ylats = activation_function.(ylats)
    g_ylats_derivative = ForwardDiff.derivative.(activation_function, ylats)
    dlogf_ylats    = - ycorr * inv(σ_ylats)
    dlogfy         = ((yobs .- μ1 - g_ylats*w1)/σ_yobs) * w1' .* g_ylats_derivative #size: (n, l1)
    gradient_ylats = dlogf_ylats + dlogfy

    return gradient_ylats  #size (n,l1)
end

# helper 2: calculate log p(z|y) to help calculate the acceptance rate
function calc_log_p_z(ylats,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)
    μ1  = weights_NN[1]
    w1 = weights_NN[2:end]
    g_ylats = activation_function.(ylats)
    logf_ylats = -0.5*sum((ycorr.^2)*inv(σ_ylats),dims=2) .- (0.5*log(prod(diag(σ_ylats))))
    logfy      = -0.5*(yobs .- μ1 - g_ylats*w1).^2 /σ_yobs .- 0.5*log(σ_yobs)
    log_p_ylats= logf_ylats + logfy

    return log_p_ylats  #size: (n,1)
end

#helper 3: one iterations of HMC to sample Z
#ycor is a temporary variable to save ycorr after reshape; ycorr is residual for latent traits
function hmc_one_iteration(nLeapfrog,ϵ,ylats_old,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)
    nobs, ntraits  = size(ylats_old)
    ylats_old = copy(ylats_old)
    ylats_new = copy(ylats_old)

    #step 1: Initiate Φ from N(0,M)
    Φ = randn(nobs, ntraits) #rand(n,Normal(0,M=1.0)), tuning parameter: M
    log_p_old = calc_log_p_z(ylats_old,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function) - 0.5*sum(Φ.^2,dims=2)  #(n,1)
    #step 2: update (ylats,Φ) from 10 leapfrog
    #2(a): update Φ
    Φ += 0.5 * ϵ * calc_gradient_z(ylats_new,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)  #(n,l1)
    for leap_i in 1:nLeapfrog
       #2(b) update latent traits
       ylats_new += ϵ * Φ  # (n,l1)
       ycorr     += ϵ * Φ  #update ycorr due to change of Z
       #(c) half step of phi
       if leap_i == nLeapfrog
           #2(c): update Φ
           Φ += 0.5 * ϵ * calc_gradient_z(ylats_new,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)
       else
           #2(a)+2(c): update Φ
           Φ += ϵ * calc_gradient_z(ylats_new,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)
       end
    end

    #Step3. acceptance rate
    log_p_new = calc_log_p_z(ylats_new,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function) - 0.5*sum(Φ.^2,dims=2) #(n,1)
    r         = exp.(log_p_new - log_p_old)  # (n,1)
    nojump    = rand(nobs) .> r  # bool (n,1)

    for i in 1:nobs
        if nojump[i]
            ylats_new[i,:] = ylats_old[i,:]
        end
    end

    return ylats_new
end


########## MT-NNMM (multiple observed traits)
#helper 1: calculate gradiant of all latent traits for all individual
function MT_calc_gradient_z(ylats,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function) #ylats=z
    μ1, w1  = weights_NN[1,:], weights_NN[2:end,:]'  #μ1:ntrais-by-1; w1: ntrait-by-nmiddle
    g_ylats = activation_function.(ylats') #ylats: nInd-by-nmiddle; g_ylats: nmiddle-by-nInd
    g_ylats_derivative = ForwardDiff.derivative.(activation_function, ylats') #g_ylats_derivative: nmiddle-by-nInd
    dlogf_ylats    = - ycorr * inv(σ_ylats) #ycorr & dlogf_ylats: nInd-by-nmiddle;
	nInd   = size(yobs,2)
    dlogfy = (w1' * inv(σ_yobs) * (yobs - μ1*ones(nInd)' - w1*g_ylats)) .* g_ylats_derivative #(nMiddle,nInd).*(nMiddle,nInd)=(nMiddle,nInd)
    gradient_ylats = dlogf_ylats + dlogfy' # (nInd,nmiddle)+(nInd,nmiddle)

    return gradient_ylats  #size: (nInd,nmiddle)
end

# helper 2: calculate log p(z|y) to help calculate the acceptance rate
function MT_calc_log_p_z(ylats,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)
    μ1, w1  = weights_NN[1,:],weights_NN[2:end,:]' #μ1:ntrais-by-1; w1: ntrait-by-nmiddle
	g_ylats = activation_function.(ylats')  #ylats: nInd-by-nmiddle; g_ylats: nmiddle-by-nInd
	nInd    = size(yobs,2)
	t = yobs - μ1*ones(nInd)' - w1*g_ylats #(ntraits,nInds)
	l = t'* σ_yobs #(nInd,ntraits)
	m = [e for e in eachrow(l)] #for each ind
	n = [e for e in eachcol(t)] #for each ind
    logf_ylats = -0.5*sum((ycorr.^2)*inv(σ_ylats),dims=2) .- 0.5*log(det(σ_ylats)) #(nInd,1)
    logfy      = -0.5*map(dot,m,n) .- 0.5*log(det(σ_yobs)) #(nInd,1)
    log_p_ylats= logf_ylats + logfy

    return log_p_ylats  #size: (nInd,1)
end

#helper 3: one iterations of HMC to sample Z
#ycor is a temporary variable to save ycorr after reshape; ycorr is residual for latent traits
# hmc_one_iteration(10,0.1, Z_missing[missingIndex,:], y[missingIndex,:]', weights_NN, I(nMiddle), σ_yobs, ycorr[missingIndex,:],activation_function)
function MT_hmc_one_iteration(nLeapfrog,ϵ,ylats_old,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)
	#ylats_old: middle nodes (nInd-by-nmiddle)
	#yobs: traits (ntraits-by-nInd)
	#weights_NN: w1 (nmiddle+1-by-ntraits)
	#ycorr: zcorr (nInd-by-nmiddle)
    nobs, ntraits  = size(ylats_old)
    ylats_old = copy(ylats_old)
    ylats_new = copy(ylats_old)

    #step 1: Initiate Φ from N(0,M)
    Φ = randn(nobs, ntraits) #rand(n,Normal(0,M=1.0)), tuning parameter: M  #1×3 Matrix Generated
    log_p_old = MT_calc_log_p_z(ylats_old,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function) - 0.5*sum(Φ.^2,dims=2)  #(n,1)

    #step 2: update (ylats,Φ) from 10 leapfrog
    #2(a): update Φ
    Φ += 0.5 * ϵ * MT_calc_gradient_z(ylats_new,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)  #(n,l1)
    for leap_i in 1:nLeapfrog
       #2(b) update latent traits
       ylats_new += ϵ * Φ  # (n,l1)
       ycorr     += ϵ * Φ  #update ycorr due to change of Z
       #(c) half step of phi
       if leap_i == nLeapfrog
           #2(c): update Φ
           Φ += 0.5 * ϵ * MT_calc_gradient_z(ylats_new,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)
       else
           #2(a)+2(c): update Φ
           Φ += ϵ * MT_calc_gradient_z(ylats_new,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function)
       end
    end

    #Step3. acceptance rate
    log_p_new = MT_calc_log_p_z(ylats_new,yobs,weights_NN,σ_ylats,σ_yobs,ycorr,activation_function) - 0.5*sum(Φ.^2,dims=2) #(n,1)
    r         = exp.(log_p_new - log_p_old)  # (n,1)
    nojump    = rand(nobs) .> r  # bool (n,1)

    for i in 1:nobs
        if nojump[i]
            ylats_new[i,:] = ylats_old[i,:]
        end
    end

    return ylats_new
end
