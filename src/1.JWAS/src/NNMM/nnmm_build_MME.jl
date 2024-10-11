"""
    build model for NMMM
"""
function nnmm_build_model(model_equations::AbstractString, 
                     ## residual variance: 
                     R = false; df = 4.0, 
                     estimate_variance=true, estimate_scale=false, 
                     constraint=false, #for multi-trait only, constraint=true means no residual covariance among traits
                     ## nnmm:
                     num_hidden_nodes = false, nonlinear_function = false, latent_traits=false, #nonlinear_function(x1,x2) = x1+x2
                     user_σ2_yobs = false, user_σ2_weightsNN = false)

    if R != false && !isposdef(map(AbstractFloat,R))
      error("The covariance matrix is not positive definite.")
    end
    if !(typeof(model_equations)<:AbstractString) || model_equations==""
      error("Model equations are wrong.\n
      To find an example, type ?build_model and press enter.\n")
    end
    if estimate_scale != false
      error("estimate scale for residual variance is not supported now.")
    end

    ############################################################################
    # Bayesian Neural Network
    ############################################################################
    if nonlinear_function != false  #NNBayes
      if latent_traits != false && length(latent_traits) != num_hidden_nodes
        error("The number of traits included in latent_traits is not $num_hidden_nodes (num_hidden_nodes)")
      end
      printstyled("Bayesian Neural Network is used with following information: \n",bold=false,color=:green)
      #NNBayes: check parameters
      num_hidden_nodes,is_fully_connected,is_activation_fcn = nnbayes_check_print_parameter(model_equations, num_hidden_nodes, nonlinear_function,latent_traits)
      #NNBayes: re-write model equations by treating hidden nodes as multiple traits
      model_equations = nnbayes_model_equation(model_equations,num_hidden_nodes,is_fully_connected)
    end
    is_nnbayes_partial = nonlinear_function != false && is_fully_connected==false #1.partial connected NN 2. fully connected NN + non-NN

    ############################################################################
    # All model terms (will be added to MME)
    ############################################################################
    #e.g., ""y2 = A+B+A*B""
    modelVec   = [strip(i) for i in split(model_equations,[';','\n'],keepempty=false)]
    nModels    = size(modelVec,1)
    if R != false && size(R,1) != nModels
      error("The residual covariance matrix is not a ",nModels," by ",nModels," matrix.")
    end
    lhsVec     = Symbol[]    #:y, phenotypes
    modelTerms = ModelTerm[] #initialization of an array of ModelTerm outside for loop
    dict       = Dict{AbstractString,ModelTerm}()
    for (m,model) = enumerate(modelVec)
      lhsRhs = split(model,"=")                  #"y2","A+B+A*B"
      lhs    = strip(lhsRhs[1])                  #"y2"
      lhsVec = [lhsVec;Symbol(lhs)]              #:y2
      rhsVec = split(strip(lhsRhs[2]),"+")       #"A","B","A*B"
      mTrms  = [ModelTerm(strip(trmStr),m,lhs) for trmStr in rhsVec]
      modelTerms  = [modelTerms;mTrms]           #a vector of ModelTerm
    end
    for trm in modelTerms          #make a dict for model terms
      dict[trm.trmStr] = trm
    end

    ############################################################################
    # Genotypes (will be added to MME)
    ############################################################################
    genotypes = []
    whichterm = 1
    for term in modelTerms
      term_symbol = Symbol(split(term.trmStr,":")[end])
      if isdefined(Main,term_symbol) #@isdefined can be used to test whether a local variable or object field is defined
        if typeof(getfield(Main,term_symbol)) == Genotypes
          term.random_type = "genotypes"
          genotypei = getfield(Main,term_symbol)
          genotypei.name = string(term_symbol)
          trait_names=[term.iTrait]
          if genotypei.name ∉ map(x->x.name, genotypes) #only save unique genotype
            genotypei.ntraits = is_nnbayes_partial ? 1 : nModels
            genotypei.trait_names = is_nnbayes_partial ? trait_names : string.(lhsVec)
            if nModels != 1
              genotypei.G.df = genotypei.G.df + nModels
            end
            if !is_nnbayes_partial && (genotypei.G.val != false || genotypei.genetic_variance.val != false)
              if size(genotypei.G.val,1) != nModels && size(genotypei.genetic_variance.val,1) != nModels
                error("The genomic covariance matrix is not a ",nModels," by ",nModels," matrix.")
              end
            end
            push!(genotypes,genotypei)
          end
        end
      end
    end

  #create mme with genotypes
  filter!(x->x.random_type != "genotypes",modelTerms) #remove "genotypes" from modelTerms
  filter!(x->x[2].random_type != "genotypes",dict)    #remove "genotypes" from dict
  

  #set scale and df for residual variance
  if nModels == 1
    scale_R = R*(df - 2)/df
    df_R    = df 
  else
    scale_R = R*(df - 1)
    df_R    = df + nModels
  end
  
  #initialize mme
  mme = MME(nModels,modelVec,modelTerms,dict,lhsVec, 
            Variance(R==false ? R : Float32.(R), #val
                     Float32(df_R),              #df
                     R==false ? R : scale_R,     #scale
                     estimate_variance, estimate_scale, constraint))
  if length(genotypes) != 0
    mme.M = genotypes #add genotypes into mme
  end

  #NNBayes:
  if nonlinear_function != false
    mme.is_fully_connected   = is_fully_connected
    mme.is_activation_fcn    = is_activation_fcn
    mme.nonlinear_function   = isa(nonlinear_function, Function) ? nonlinear_function : nnbayes_activation(nonlinear_function)
    mme.latent_traits        = latent_traits
    if user_σ2_yobs != false && user_σ2_weightsNN != false
      mme.σ2_yobs         = user_σ2_yobs      #variance of observed phenotype σ2_yobs is fixed as user_σ2_yobs
      mme.σ2_weightsNN    = user_σ2_weightsNN #variance of neural network weights between omics and phenotype σ2_weightsNN is fixed as user_σ2_weightsNN
      mme.fixed_σ2_NN     = true
      printstyled(" - Variances of phenotype and neural network weights are fixed.\n",bold=false,color=:green)
    end
  end

  return mme
end