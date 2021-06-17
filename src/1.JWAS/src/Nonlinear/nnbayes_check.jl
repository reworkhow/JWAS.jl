#Below function is to check parameters for NNBayes and print information
function nnbayes_check_print_parameter(num_latent_traits,nonlinear_function,activation_function)
    printstyled("Bayesian Neural Network is used with follwing information: \n",bold=false,color=:green)

    #part1: fully/partial-connected NN
    if typeof(num_latent_traits) == Int64         #fully-connected. e.g, num_latent_traits=5
        printstyled(" - Neural network:         fully connected neural network \n",bold=false,color=:green)
        printstyled(" - Number of hidden nodes: $num_latent_traits \n",bold=false,color=:green)
    elseif num_latent_traits == false   #partial-connected.
        printstyled(" - Neural network:         partially connected neural network \n",bold=false,color=:green)
    else
        error("Please check you number of latent traits")
    end

    #part2: activation function/user-defined non-linear function
    if nonlinear_function == "Neural Network" #NN
        if activation_function in ["tanh","sigmoid","relu","leakyrelu","linear"]
            printstyled(" - Activation function:    $activation_function.\n",bold=false,color=:green)
            printstyled(" - Sampler:                Hamiltonian Monte Carlo. \n",bold=false,color=:green)
        else
            error("Please select the activation function from tanh/sigmoid/relu/leakyrelu/linear")
        end
    elseif isa(nonlinear_function, Function) #user-defined nonlinear function. e.g, CropGrowthModel()
        if activation_function == false
            printstyled(" - Nonlinear function:     user-defined nonlinear_function for the relationship between hidden nodes and observed trait is used.\n",bold=false,color=:green)
            printstyled(" - Sampler:                Matropolis-Hastings.\n",bold=false,color=:green)
        else
            error("activation function is not allowed for user-defined nonlinear function")
        end
    else
        error("nonlinear_function can only be Neural Network or a user-defined nonlinear function")
    end
end


#Below function is to re-phase modelm for NNBayes
function nnbayes_model_equation(model_equations,num_latent_traits)

    lhs, rhs = strip.(split(model_equations,"="))
    model_equations = ""

    if typeof(num_latent_traits) == Int64   #fully-connected
      # old: y=intercept+geno
      # new: y1=intercept+geno;y2=intercept+geno
      for i = 1:num_latent_traits
        model_equations = model_equations*lhs*string(i)*"="*rhs*";"
      end
    elseif num_latent_traits == false      #partially-connected
      # old: y=intercept+geno1+geno2
      # new: y1= intercept+geno1;y2=intercept+geno2
      rhs_split=strip.(split(rhs,"+"))
      geno_term=[]
      for i in rhs_split
        if isdefined(Main,Symbol(i)) && typeof(getfield(Main,Symbol(i))) == Genotypes
            push!(geno_term,i)
        end
      end
      non_gene_term = filter(x->x ∉ geno_term,rhs_split)
      non_gene_term = join(non_gene_term,"+")

      for i = 1:length(geno_term)
        model_equations = model_equations*lhs*string(i)*"="*non_gene_term*"+"*geno_term[i]*";"
      end
    end
    model_equations = model_equations[1:(end-1)]
end


# below function is to check whether the loaded genotype matches the model equation
function nnbayes_check_nhiddennode(num_latent_traits,mme)
    if typeof(num_latent_traits) == Int64 #fully-connected. e.g, num_latent_traits=5
        if length(mme.M)>1
            error("fully-connected NN only allow one genotype; num_latent_traits is not allowed in partial-connected NN ")
        end
    elseif num_latent_traits == false   #partial-connected.
        if length(mme.M)==1
            error("partial-connected NN requirs >1 genotype group")
        else
            num_latent_traits = length(mme.M)
            printstyled(" - Number of hidden nodes: $num_latent_traits \n",bold=false,color=:green)
        end #Note, if only geno1 & geno2 are loaded by get_genotypes, but there is "geno3" in equation, then geno3 will be treated like age.
    end
end



# below function is to define the activation function for neural network
function nnbayes_activation(activation_function)
      if activation_function == "tanh"
         mytanh(x) = tanh(x)
         return mytanh
      elseif activation_function == "sigmoid"
         mysigmoid(x) = 1/(1+exp(-x))
         return mysigmoid
      elseif activation_function == "relu"
         myrelu(x) = max(0, x)
         return myrelu
      elseif activation_function == "leakyrelu"
         myleakyrelu(x) = max(0.01x, x)
         return myleakyrelu
      elseif activation_function == "linear"
         mylinear(x) = x
         return mylinear
      else
          error("invalid actication function")
      end
end


# below function is to modify mme from multi-trait model to multiple single trait models
# coded by Hao
function nnbayes_mega_trait(mme)
    #mega_trait
    if mme.nModels == 1
        error("more than 1 trait is required for MegaLMM analysis.")
    end
    mme.MCMCinfo.constraint = true

    ##sample from scale-inv-⁠χ2, not InverseWishart
    mme.df.residual  = mme.df.residual - mme.nModels
    mme.scaleR       = diag(mme.scaleR/(mme.df.residual - 1))*(mme.df.residual-2)/mme.df.residual #diag(R_prior_mean)*(ν-2)/ν
    if mme.M != 0
        for Mi in mme.M
            Mi.df        = Mi.df - mme.nModels
            Mi.scale    = diag(Mi.scale/(Mi.df - 1))*(Mi.df-2)/Mi.df
        end
    end

end


# below function is to modify essential parameters for partial connected NN
function nnbayes_partial_para_modify(mme)
    mme.nnbayes_partial=true
    for Mi in mme.M
      Mi.ntraits=1
    end
end


# below function is to modify essential parameters for partial connected NN
function nnbayes_partial_para_modify2(mme)
    for Mi in mme.M
      Mi.scale = Mi.scale[1]
      Mi.G = Mi.G[1,1]
      Mi.genetic_variance=Mi.genetic_variance[1,1]
    end
end


# below function is to modify essential parameters for partial connected NN
function nnbayes_partial_para_modify3(mme)
    for Mi in mme.M
      Mi.meanVara  = Mi.meanVara[1]
      Mi.meanVara2 = Mi.meanVara2[1]
      Mi.meanScaleVara = Mi.meanScaleVara[1]
      Mi.meanScaleVara2 = Mi.meanScaleVara2[1]
      Mi.π = Mi.π[1]
      Mi.mean_pi = Mi.mean_pi[1]
      Mi.mean_pi2 = Mi.mean_pi2[1]
    end
end
