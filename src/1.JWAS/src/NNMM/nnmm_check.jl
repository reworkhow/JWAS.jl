#Below function is to check parameters for NNBayes and print information
function nnbayes_check_print_parameter(model_equations,num_hidden_nodes,nonlinear_function,latent_traits)
    ## Step1. determine the neural network architecture
    # (i) count the number of genotypes
    model_terms = Symbol.(strip.(split(split(model_equations,"=")[2],"+")))
    ngeno = 0
    for i in model_terms
        if isdefined(Main,i) && typeof(getfield(Main,i)) == Genotypes
            ngeno+=1
        end
    end
    if ngeno==0
        error("please load genotypes.")
    end

    # (ii) determine whether partial/fully connected neural network
    if isa(nonlinear_function, Function)   #e.g., nonlinear_function is PGM
        nargs_nonlinear = first(methods(nonlinear_function)).nargs-1 #e.g., PGM(x1,x2) has 2 arguments
        is_activation_fcn = false
        if ngeno == 1    # fully-connected (equivalent to y=geno1+geno1+geno1 (memory efficient))
            is_fully_connected = true
            num_hidden_nodes = nargs_nonlinear
        elseif ngeno == nargs_nonlinear # partial-connected
            is_fully_connected = false
            num_hidden_nodes = ngeno
        else
            error(" number of arguments in nonlinear_function ≠ number of loaded genotypes.")
        end
    elseif nonlinear_function in ["tanh","sigmoid","relu","leakyrelu","linear"]
        is_activation_fcn = true
        if num_hidden_nodes != false && typeof(num_hidden_nodes) == Int64
            num_hidden_nodes = num_hidden_nodes
            if ngeno == 1 # fully-connected
                is_fully_connected = true
            elseif ngeno  == num_hidden_nodes #partial-connected
                is_fully_connected = false
            else
                error("number of loaded genotypes ≠ number of hidden nodes")
            end
        elseif num_hidden_nodes == false  #partial
            is_fully_connected = false
            num_hidden_nodes   = ngeno
        else
            error("The num_hidden_nodes should be an interger.")
        end
    else
        error("invalid nonlinear_function.")
    end

    # (iii) print NNBayes info
    if is_fully_connected == true
        printstyled(" - Neural network:         fully connected neural network. \n",bold=false,color=:green)
    elseif is_fully_connected == false
        printstyled(" - Neural network:         partially connected neural network. \n",bold=false,color=:green)
    else
        error("error")
    end
    printstyled(" - Number of hidden nodes: $num_hidden_nodes. \n",bold=false,color=:green)

    if is_activation_fcn==true  #NN with activation function
        printstyled(" - Nonlinear function:     $nonlinear_function.\n",bold=false,color=:green)
        printstyled(" - Sampler:                Hamiltonian Monte Carlo. \n",bold=false,color=:green)
    elseif is_activation_fcn==false #user-defined nonlinear function. e.g, CropGrowthModel()
        printstyled(" - Nonlinear function:     user-defined nonlinear_function for the relationship between hidden nodes and observed trait is used.\n",bold=false,color=:green)
        printstyled(" - Sampler:                Matropolis-Hastings.\n",bold=false,color=:green)
    else
        error("invalid nonlinear_function")
    end

    if latent_traits!=false
        nomics=length(latent_traits)
        printstyled(" - Omics data:             $nomics intermediate omics features are included.\n",bold=false,color=:green)
    end

    return num_hidden_nodes,is_fully_connected,is_activation_fcn
end


#Below function is to re-phase modelm for NNBayes
function nnbayes_model_equation(model_equations,num_hidden_nodes,is_fully_connected)

    lhs, rhs = strip.(split(model_equations,"="))
    model_equations = ""

    if is_fully_connected == true   #fully-connected
      # old: y=intercept+geno
      # new: y1=intercept+geno;y2=intercept+geno
      for i = 1:num_hidden_nodes
        model_equations = model_equations*lhs*string(i)*"="*rhs*";"
      end
    elseif is_fully_connected == false   #partially-connected
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
    model_equations = model_equations[1:(end-1)] #remove the semicolon at the end
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
    mme.R.df    = mme.R.df - mme.nModels
    mme.R.scale = diag(mme.R.scale/(mme.R.df - 1))*(mme.R.df-2)/mme.R.df #diag(R_prior_mean)*(ν-2)/ν
    if mme.M != 0
        for Mi in mme.M
            Mi.G.df    = Mi.G.df - mme.nModels
            Mi.G.scale = diag(Mi.G.scale/(Mi.G.df - 1))*(Mi.G.df-2)/Mi.G.df
        end
    end

end

# below function is to modify essential parameters for partial connected NN
function nnbayes_partial_para_modify2(mme)
    for Mi in mme.M
      Mi.G.scale = Mi.G.scale[1]
      Mi.G.val = Mi.G.val[1,1]
      Mi.genetic_variance.val=Mi.genetic_variance.val[1,1]
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



#below function is to initialize missing omics with yobs
function nnlmm_initialize_missing_with_mean(mme,df)
    #define mme.yobs here because ungenotyped inds has been removed from df
    mme.yobs = df[!,mme.yobs_name]  # mme.yobs = DataFrames.recode(df[!,mme.yobs_name], missing => 0.0)  #e.g., mme.lhsVec=[:y1,:y2]
    # if mme.latent_traits != false  #NN-LMM-Omics
      #save omics data missing pattern
      mme.missingPattern = .!ismissing.(Matrix(df[!,mme.lhsVec]))
      #replace missing data with values in yobs
      for i in mme.lhsVec      #for each omics feature
        all_missing_omics_i = all(ismissing.(df[!,i])) #whether all omics i are missing
        if !all_missing_omics_i #not all omics i are missing
            mean_omics_i = mean(filter(isfinite,skipmissing(df[!,i]))) #mean of non-missing omics i
        else #all omics i are missing, change col type from Missing to union of Missing and Float64
            df[!,i] = convert(Vector{Union{Missing,Float64}}, df[!,i])
        end
        for j in 1:size(df,1)  #for each observation
          if ismissing(df[j,i])
            if !all_missing_omics_i
                df[j,i]=mean_omics_i
            else #all omics i are missing
                #if yobs[j] is missing, remove the row from df
                if ismissing(mme.yobs[j])
                    error("this individual has no omics data and no phenotype data: $(df[j,:ID]), please remove this individual.")
                else
                    df[j,i]=mme.yobs[j]
                end
            end
          end
        end
      end
    # else  #NN-Bayes with hidden nodes (G3 paper)
    #   #all omics should be missing, the missingPattern should be all 0
    #   #but we already set y1,...,y5 as yobs, so we have to build missingPattern
    #   # byhand.
    #   mme.missingPattern = .!ismissing.(Array{Missing}(missing, size(df[!,mme.lhsVec])))
    # end
    # add indicators for individuals with full omics data, so their omics won't be sampled
    n_observed_omics = sum(mme.missingPattern,dims=2) #number of observed omics for each ind
    n_omics          = length(mme.lhsVec)             #number of omics
    full_omics       = n_observed_omics .== n_omics   #indicator for ind with full omics
    mme.incomplete_omics    = vec(.!full_omics)              #indicator for ind with no/partial omics
end


function nnmm_print_info_input_to_middle_layer(mme)
        nThread = Threads.nthreads()
        is_mega_trait = mme.R.constraint==true && mme.M[1].G.constraint==true #is_mega_trait when no residual and marker effect covariances 
        if is_mega_trait
            printstyled(" - Bayesian Alphabet:                multiple independent single-trait Bayesian models are used to sample marker effect. \n",bold=false,color=:green)
            printstyled(" - Multi-threaded parallelism:       $nThread threads are used to run single-trait models in parallel. \n",bold=false,color=:green)
        else
            printstyled(" - Bayesian Alphabet:                multi-trait Bayesian models are used to sample marker effect. \n",bold=false,color=:green)
        end
end