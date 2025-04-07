function runNNMM(layers, equations;
                #MCMC
                chain_length::Integer             = 100,
                burnin::Integer                   = 0,
                output_samples_frequency::Integer = chain_length>1000 ? div(chain_length,1000) : 1,
                update_priors_frequency::Integer  = 0,
                #Genomic Prediction
                outputEBV                       = true,
                output_heritability             = true,
                prediction_equation             = false,
                #MISC
                seed                            = false,
                printout_model_info             = true,
                printout_frequency              = chain_length+1,
                big_memory                      = false,
                double_precision                = false,
                #MCMC samples (defaut to marker effects and hyperparametes (variance componets))
                output_folder                     = "nnmm_results",
                output_samples_for_all_parameters = false) 

    ############################################################################
    # Step1. read genotypes
    #  note: need to read genotypes here because "method" is defined in equations
    ############################################################################
    printstyled("Starte reading genotypes.......\n",bold=false,color=:green)
    #arguments for data
    file = layers[1].data_path
    separator = layers[1].separator
    header = layers[1].header
    quality_control = layers[1].quality_control
    MAF=layers[1].MAF
    missing_value=layers[1].missing_value
    center = layers[1].center
    
    #arguments for equations
    G    = equations[1].G
    method = equations[1].method
    Pi = equations[1].Pi
    estimatePi = equations[1].estimatePi
    G_is_marker_variance = equations[1].G_is_marker_variance
    df = equations[1].df
    estimate_variance = equations[1].estimate_variance
    estimate_scale = equations[1].estimate_scale
    constraint = equations[1].constraint
    if constraint != true
        error("please set constraint=true from 1->2 layer for multiple single-trait analysis.")
    end
    starting_value = equations[1].starting_value

    geno = get_genotypes(file, G;
                        ## method:
                        method = method, Pi = Pi, estimatePi = estimatePi, 
                        ## variance:
                        G_is_marker_variance = G_is_marker_variance, df = df,
                        estimate_variance = estimate_variance, estimate_scale = estimate_scale,
                        constraint = constraint, #for multi-trait only, constraint=true means no genetic covariance among traits
                        ## format:
                        separator = separator, header = header, double_precision = false,
                        ## quality control:
                        quality_control = quality_control, MAF = MAF, missing_value = missing_value,
                        ## others:
                        center = center, starting_value = starting_value)
    #save omics data in 2nd layer
    layers[1].data=geno

    printstyled("End reading genotypes.......\n",bold=false,color=:green)
    
    ############################################################################
    # Step2. read omics
    ############################################################################
    printstyled("Start reading omics (please replace `marker` with `omics` before finishing reading omics).......\n",bold=false,color=:green)
    printstyled(" Note that `center`, `quality_control`, and `MAF` are set to false here forcefully, so please use pre-processed omics data \n",bold=false,color=:green)
    #arguments for data
    file            = layers[2].data_path
    separator       = layers[2].separator
    header          = layers[2].header
    quality_control = layers[2].quality_control = false #QC=false -> missing omics will not be replaced by its mean
    MAF             = layers[2].MAF             = false
    center          = layers[2].center          = false
    missing_value   = layers[2].missing_value

    #arguments for equations
    G = equations[2].G
    method = equations[2].method
    if method == "GBLUP"
        error("GBLUP is not supported from 2->3 layers. Try BayesC, RR-BLUP, etc..")
    end
    Pi                   = equations[2].Pi
    estimatePi           = equations[2].estimatePi
    G_is_marker_variance = equations[2].G_is_marker_variance
    df                   = equations[2].df
    estimate_variance    = equations[2].estimate_variance
    estimate_scale       = equations[2].estimate_scale
    constraint           = equations[2].constraint = false
    printstyled(" Note that `center`, `quality_control`, and `MAF` are set to false here forcefully, so please use pre-processed omics data \n",bold=false,color=:green)
    starting_value       = equations[2].starting_value

    omcis = get_genotypes(file, G;
                    ## method:
                    method = method, Pi = Pi, estimatePi = estimatePi, 
                    ## variance:
                    G_is_marker_variance = G_is_marker_variance, df = df,
                    estimate_variance = estimate_variance, estimate_scale = estimate_scale,
                    constraint = false, #for multi-trait only, constraint=true means no genetic covariance among traits
                    ## format:
                    separator = separator, header = header, double_precision = false,
                    ## quality control:
                    quality_control = false, MAF = false, missing_value = missing_value,
                    ## others:
                    center = false, starting_value = starting_value)
    #reset some parameters for omics
    omics.alleleFreq = false
    omics.sum2pq     = false
    #save omics data in 2nd layer
    layers[2].data=omcis
    
    printstyled("End reading omics.......\n",bold=false,color=:green)
    
    ############################################################################
    #Step3. check input equations
    ############################################################################
    printstyled("Start check input data. \n",bold=false,color=:green)
    #check the number of equations
    if length(equations) != 2
        error("The number of equations must be 2: (1) first equation from 1st layer to 2nd layer; (2) second equation from 2nd layer to 3rd layer.")
    end
    #check the order of equations
    if equations[1].to_layer_name != equations[2].from_layer_name
        error("The order of equations is wrong: (1) first equation from 1st layer to 2nd layer; (2) second equation from 2nd layer to 3rd layer.")
    end
    #check the order of layers, and layers names
    if layers[1].layer_name != equations[1].from_layer_name || layers[2].layer_name != equations[1].to_layer_name || layers[2].layer_name != equations[2].from_layer_name || layers[3].layer_name != equations[2].to_layer_name
        error("The order of layers is wrong, or the layer name is not consistent with the equation.")
    end
    #make sure the layer name is used in equations, and only appears once
    equ_left_tmp, equ_right_tmp = strip.(split(equations[1].equation,"=")) # "omics = intercept + geno"
    n_occurrences = length(split(equ_right_tmp, layers[1].layer_name)) - 1 # how many times "geno" appears in the righ of equation, must=1
    if equ_left_tmp != layers[2].layer_name || n_occurrences != 1
        error("The layer name in equation is not consistent with the layer name, and please only include it once in your equation.")
    end
    equ_left_tmp, equ_right_tmp = strip.(split(equations[2].equation,"=")) # "phenotypes = intercept + omics"
    n_occurrences = length(split(equ_right_tmp, layers[2].layer_name)) - 1 # how many times "omics" appears in the right of equation, must=1
    if equ_left_tmp != layers[3].layer_name || n_occurrences != 1
        error("The layer name in equation is not consistent with the layer name, and please only include it once in your equation.")
    end

    #check partial_connect_structure
    is_fully_connected_12 = true #_12 means from 1st layer to 2nd layer
    if equations[1].partial_connect_structure != false
        #make a marker effect model is used for partial connected 1->2 layer
        if equations[1].method=="GBLUP"
            error("The partial_connect_structure is not allowed for GBLUP model.")
        end
        #check if all elements are either 0 or 1, change to sparse matrix
        if all(equations[1].partial_connect_structure .== 0 .|| equations[1].partial_connect_structure .== 1)
            #change to sparse boolean matrix
            equations[1].partial_connect_structure = sparse(equations[1].partial_connect_structure .== 1) #change to sparse Boolean matrix
        else
            error("Elements in partial_connect_structure should be either 0 or 1.")
        end
        #check the dimention
        # - number of snps in the 1st layer = layers[1].data.nMarkers
        # - number of omics in the 2nd layer = layers[2].data.nMarkers
        if (layers[1].data.nMarkers,layers[2].data.nMarkers) != size(equations[1].partial_connect_structure)
            error("From 1->2, the dimention of partial_connect_structure is wrong. Correct dimention should be : #nodes_layer1 - by - #nodes_layer2.")
        end
        is_fully_connected_12 = false
    end

    is_fully_connected_23 = true #_23 means from 1st layer to 2nd layer
    if equations[2].partial_connect_structure != false
        error("The partial_connect_structure is currently only allowed from 1->2 layers.")
        is_fully_connected_23 = false
    end
    #check actication function
    if isa(equations[2].activation_function, Function)  #e.g., user-defined activation_function=Pig_Growth_Model(o1, o2, o3)
        num_args = first(methods(equations[2].activation_function)).nargs-1 #number of arguments in activation_function, should = number of nodes in the 2nd layer
        if num_args != layers[2].data.nMarkers #number of omics in the 2nd layer = layers[2].data.nMarkers
            error("The number of arguments in user-defined activation_function ≠ number of nodes in the 2nd layer.")
        end
        is_activation_fcn = false
    elseif equations[2].activation_function in ["tanh","sigmoid","relu","leakyrelu","linear"]
        is_activation_fcn = true
    else
        error("invalid activation_function.")
    end
    #check whether the 1st layer is genotype
    if typeof(layers[1].data) != Genotypes
        error("currently, the 1st layer should be genotypes in NNMM.")
    end

    ############################################################################
    #Step4. print NNMM info
    ############################################################################
    printstyled("******************** NNMM information ******************** \n",bold=false,color=:green)
    if is_fully_connected_12 == true
        printstyled(" - From 1st->2nd layers:             fully connected network. \n",bold=false,color=:green)
    elseif is_fully_connected_12 == false
        printstyled(" - From 1st->2nd layers:             partially connected network. \n",bold=false,color=:green)
    else
        error("error")
    end
        printstyled(" - Number of nodes in 2nd layer:     $(layers[2].data.nMarkers). \n",bold=false,color=:green)
    #check if data in 2nd layer are all missing
    if sum(ismissing.(layers[2].data.data)) == size(layers[2].data.data,1) * size(layers[2].data.data,2)
        printstyled("  - no omics data!! 2nd layer are all hidden nodes.\n",bold=false,color=:green)
    end

    if is_activation_fcn==true  #NN with activation function such as tanh, sigmoid
        printstyled(" - Activation function in 2nd layer: $(equations[2].activation_function).\n",bold=false,color=:green)
        printstyled(" - Sample missing data in 2nd layer: Hamiltonian Monte Carlo. \n",bold=false,color=:green)
    elseif is_activation_fcn==false #user-defined function. e.g, Pig_Growth_Model(o1, o2, o3)
        printstyled(" - Nonlinear function in 2nd layer:  user-defined nonlinear_function for the relationship between middle layer and output layer is used.\n",bold=false,color=:green)
        printstyled(" - Sample missing data in 2nd layer: Matropolis-Hastings.\n",bold=false,color=:green)
    else
        error("invalid nonlinear_function")
    end

    
    ############################################################################
    # Step5. build model equation (1->2): replace "omics" with omics feature names
    #  - old: omics=intercept+geno
    #  - new: o1=intercept+geno;o2=intercept+geno
    ############################################################################
    equation_12 = ""
    _, equ_right = strip.(split(equations[1].equation,"=")) # "intercept + x1 + geno" 
    for i in layers[2].data.featureID
        equation_12 = equation_12 * i * "=" * equ_right * ";"
    end
    equation_12 = equation_12[1:(end-1)] #remove the ";" at the end
    equations[1].equation = equation_12 #update the equation
    ############################################################################
    # Step6. build model equation (2->3): replace "phenotypes" with trait names
    #  - old: phenotypes=intercept+omics
    #  - new: single-trait: y1=intercept+omics
    #         multi-trait:  y1=intercept+omics;y2=intercept+omics
    ############################################################################
    equation_23 = ""
    _, equ_right = strip.(split(equations[2].equation,"=")) # "intercept + omics"
    for i in layers[3].data.featureID
        equation_23 = equation_23 * i * "=" * equ_right * ";"
    end
    equation_23 = equation_23[1:(end-1)] #remove the ";" at the end
    equations[2].equation = equation_23 #update the equation

    mme_all = []
    ############################################################################
    # 1->2: All model terms (will be added to MME)
    # same as build_model()
    ############################################################################
    model_equations=equations[1].equation
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
    # 1->2: Genotypes (will be added to MME)
    ############################################################################
    genotypes = []
    whichterm = 1
    for term in modelTerms
        term_symbol = Symbol(split(term.trmStr,":")[end])
        println(term_symbol)
        if term_symbol == Symbol(layers[1].layer_name) #genotypes
          term.random_type = "genotypes"
          genotypei = layers[1].data
          genotypei.name = string(term_symbol)
          trait_names=[term.iTrait]
          if genotypei.name ∉ map(x->x.name, genotypes) #only save unique genotype
            genotypei.ntraits = nModels
            genotypei.trait_names = string.(lhsVec)
            if nModels != 1
              genotypei.G.df = genotypei.G.df + nModels
            end
            if (genotypei.G.val != false || genotypei.genetic_variance.val != false)
              if size(genotypei.G.val,1) != nModels && size(genotypei.genetic_variance.val,1) != nModels
                error("The genomic covariance matrix is not a ",nModels," by ",nModels," matrix.")
              end
            end
            push!(genotypes,genotypei)
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
    
    mme.is_fully_connected   = is_fully_connected_12
    
    append!(mme_all,mme)

    ############################################################################
    # 2->3: All model terms (will be added to MME)
    # same as build_model()
    ############################################################################
    model_equations=equations[2].equation
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
    # 2->3: Omics (will be added to MME)
    ############################################################################
    genotypes = []
    whichterm = 1
    for term in modelTerms
        term_symbol = Symbol(split(term.trmStr,":")[end])
        if term_symbol in layers[2].data.featureID
          term.random_type = "omics"
          genotypei = layers[2].data
          genotypei.name = string(term_symbol)
          trait_names=[term.iTrait]
          if genotypei.name ∉ map(x->x.name, genotypes) #only save unique genotype
            genotypei.ntraits = nModels
            genotypei.trait_names = string.(lhsVec)
            if nModels != 1
              genotypei.G.df = genotypei.G.df + nModels
            end
            if (genotypei.G.val != false || genotypei.genetic_variance.val != false)
              if size(genotypei.G.val,1) != nModels && size(genotypei.genetic_variance.val,1) != nModels
                error("The genomic covariance matrix is not a ",nModels," by ",nModels," matrix.")
              end
            end
            push!(genotypes,genotypei)
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
    
    mme.is_fully_connected   = is_fully_connected_12
    
    append!(mme_all,mme)





    ############################################################################
    # Neural Network
    ############################################################################
    is_nnbayes_partial = (mme.nonlinear_function != false && mme.is_fully_connected==false)
    if mme.nonlinear_function != false #modify data to add phenotypes for hidden nodes
        mme.yobs_name=Symbol(mme.lhsVec[1]) #e.g., lhsVec=[:y1,:y2,:y3], a number label has been added to original trait name in nnbayes_model_equation(),
        yobs = df[!,Symbol(string(mme.yobs_name)[1:(end-1)])]  # e.g., change :y1 -> :y
        for i in mme.lhsVec  #e.g., lhsVec=[:y1,:y2,:y3]
            df[!,i]= yobs
        end
        ######################################################################
        #mme.lhsVec and mme.M[1].trait_names default to empirical trait name
        #with prefix 1, 2... , e.g., height1, height2...
        #if data for latent traits are included in the dataset, column names
        #will be used as shown below.e.g.,
        #mme.latent_traits=["gene1","gene2"],  mme.lhsVec=[:gene1,:gene2] where
        #"gene1" and "gene2" are columns in the dataset.
        ######################################################################
        if mme.latent_traits != false
            #change lhsVec to omics gene name
            mme.lhsVec = Symbol.(mme.latent_traits) # [:gene1, :gene2, ...]
            #rename genotype names
            mme.M[1].trait_names=mme.latent_traits
            #change model terms for partial-connected NN
            if is_nnbayes_partial
                for i in 1:mme.nModels
                    mme.M[i].trait_names=[mme.latent_traits[i]]
                end
            end
        end
        mme.R.constraint=true
        for Mi in mme.M
            Mi.G.constraint=true
        end
        nnmm_print_info_input_to_middle_layer(mme)
    end
    ############################################################################
    # Set a seed in the random number generator
    ############################################################################
    if seed != false
        Random.seed!(seed)
    end
    #when using multi-thread, make sure the results are reproducible for users
    nThread = Threads.nthreads()
    if nThread>1
        Threads.@threads for i = 1:nThread
              Random.seed!(seed+i)
        end
    end
    ############################################################################
    # Create an folder to save outputs
    ############################################################################
    myfolder,folderi = output_folder, 1
    while ispath(output_folder)
        printstyled("The folder $output_folder already exists.\n",bold=false,color=:red)
        output_folder = myfolder*string(folderi)
        folderi += 1
    end
    mkdir(output_folder)
    printstyled("The folder $output_folder is created to save results.\n",bold=false,color=:green)
    ############################################################################
    # Save MCMC argumenets in MCMCinfo
    ############################################################################
    mme.MCMCinfo = MCMCinfo(heterogeneous_residuals,
                   chain_length,burnin,output_samples_frequency,
                   printout_model_info,printout_frequency, single_step_analysis,
                   fitting_J_vector,missing_phenotypes,
                   update_priors_frequency,outputEBV,output_heritability,prediction_equation,
                   seed,double_precision,output_folder,false)
    ############################################################################
    # Check 1)Arguments; 2)Input Pedigree,Genotype,Phenotypes,
    #       3)output individual IDs; 4)Priors  5)Prediction equation
    #       in the Whole Dataset.
    ############################################################################
    errors_args(mme)       #check errors in function arguments
    df=check_pedigree_genotypes_phenotypes(mme,df,pedigree)
    if mme.nonlinear_function != false #NN-LMM
        #initiliza missing omics data  (after check_pedigree_genotypes_phenotypes() because non-genotyped inds are removed)
        nnlmm_initialize_missing(mme,df)
    end
    prediction_setup(mme)  #set prediction equation, defaulting to genetic values
    check_outputID(mme)    #check individual of interest for prediction
    df_whole,train_index = make_dataframes(df,mme)
    set_default_priors_for_variance_components(mme,df_whole)  #check priors (set default priors)

    if mme.M!=0
        set_marker_hyperparameters_variances_and_pi(mme)
    end
    ############################################################################
    # Adhoc functions
    ############################################################################
    #save MCMC samples for all parameters (?seperate function user call)
    if output_samples_for_all_parameters == true
        allparameters=[term[2] for term in split.(keys(mme.modelTermDict),":")]
        allparameters=unique(allparameters)
        for parameter in allparameters
            outputMCMCsamples(mme,parameter)
        end
    end
    # Double Precision
    if double_precision == true
        if mme.M != 0
            for Mi in mme.M
                Mi.genotypes = map(Float64,Mi.genotypes)
                Mi.G.val         = map(Float64,Mi.G.val)
                Mi.α         = map(Float64,Mi.α)
            end
        end
        for random_term in mme.rndTrmVec
            random_term.Vinv  = map(Float64,random_term.Vinv)
            random_term.GiOld.val = map(Float64,random_term.GiOld.val)
            random_term.GiNew.val = map(Float64,random_term.GiNew.val)
            random_term.Gi.val    = map(Float64,random_term.Gi.val)
        end
        mme.Gi = map(Float64,mme.Gi)
    end

    #constraint on covariance matrix
    if mme.R.constraint==true
        R_constraint!(mme) #modify mme.R.df and mme.R.scale; scale is a diagonal matrix
    end
    if mme.M != 0 && mme.M[1].G.constraint==true
        G_constraint!(mme) #modify Mi.G.df and Mi.G.scale; scale is a diagonal matrix
    end

    # NNBayes: modify parameters for partial connected NN
    if is_nnbayes_partial
        nnbayes_partial_para_modify2(mme)
    end
    ############################################################################
    #Make incidence matrices and genotype covariates for training observations
    #and individuals of interest
    ############################################################################
    #make incidence matrices (non-genomic effects) (after SSBRrun for ϵ & J)
    df=make_incidence_matrices(mme,df_whole,train_index)
    #align genotypes with 1) phenotypes IDs; 2) output IDs.
    if mme.M != false
        align_genotypes(mme,output_heritability,single_step_analysis)
    end
    # initiate Mixed Model Equations and check starting values
    init_mixed_model_equations(mme,df,starting_value)
    ############################################################################
    # MCMC
    ############################################################################
    describe(mme)
    mme.output=nnmm_MCMC_BayesianAlphabet(mme,df)

    ############################################################################
    # Save output to text files
    ############################################################################
    for (key,value) in mme.output
      CSV.write(output_folder*"/"*replace(key," "=>"_")*".txt",value)
    end

    if mme.M != 0
        for Mi in mme.M
            if Mi.name == "GBLUP"
                mv("marker_effects_variance_"*Mi.name*".txt","genetic_variance(REML)_"*Mi.name*".txt")
            end
        end
    end
    printstyled("\n\nThe version of Julia and Platform in use:\n\n",bold=true)
    versioninfo()
    printstyled("\n\nThe analysis has finished. Results are saved in the returned ",bold=true)
    printstyled("variable and text files. MCMC samples are saved in text files.\n\n\n",bold=true)

    return mme.output

end