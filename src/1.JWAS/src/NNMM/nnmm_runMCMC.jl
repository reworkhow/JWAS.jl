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
    # Step1. read genotypes in layer 1
    #  note: need to read genotypes here because "method" is defined in equations
    ############################################################################
    if layers[1].data != [] #not empty, user should re-run layers function to reset data stored in layer 1
        error("There is data already stored in layer 1, please re-run the code to reset layers.")
    end

    printstyled("********************************************************** \n",bold=false,color=:green)
    printstyled("Start reading genotypes.......\n",bold=false,color=:green)
    #arguments for data
    file_paths = layers[1].data_path #file is a string for full connected; or a vector of strings for partial connected   
    #full connected, make it a vector of 1 element
    if typeof(file_paths) == String
        file_paths=[file_paths]
    end
    separator = layers[1].separator
    header = layers[1].header
    quality_control = layers[1].quality_control
    MAF=layers[1].MAF
    missing_value=layers[1].missing_value
    center = layers[1].center
    
    #arguments for equations for 1->2 layer
    G    = equations[1].G
    method = equations[1].method
    Pi = equations[1].Pi
    estimatePi = equations[1].estimatePi
    
    G_is_marker_variance = equations[1].G_is_marker_variance
    df_G = equations[1].df_G
    estimate_variance_G = equations[1].estimate_variance_G
    estimate_scale_G = equations[1].estimate_scale_G
    constraint_G = equations[1].constraint_G #by default is true
    if constraint_G != true
        error("please set constraint=true from 1->2 layer for multiple single-trait analysis.")
    end
    starting_value = equations[1].starting_value
    if starting_value != false
        error("user-defined starting value is not supported from 1->2 layer.")
    end
    
    #partial connected if multiple genotype files are provided
    is_fully_connected_12 = true #_12 means from 1st layer to 2nd layer
    is_nnbayes_partial = false
    if length(file_paths) > 1
        equations[1].partial_connect_structure = true
        is_fully_connected_12 = false
        is_nnbayes_partial = true
        printstyled("Partial connected neural network is used from 1->2 layer for multiple genotype files.\n",bold=false,color=:green)
    end
    
    genotype_name_list = []
    for i in 1:length(file_paths)
        file_path = file_paths[i]
        println()
        printstyled("Reading $(i)th genotypes from $file_path \n",bold=false,color=:green)
        # geno is of type "genotypes"
        geno = nnmm_get_genotypes(file_path, G;
                            ## method:
                            method = method, Pi = Pi, estimatePi = estimatePi, 
                            ## variance:
                            G_is_marker_variance = G_is_marker_variance, df = df_G,
                            estimate_variance = estimate_variance_G, estimate_scale = estimate_scale_G,
                            constraint = constraint_G, #for multi-trait only, constraint=true means no genetic covariance among traits
                            ## format:
                            separator = separator, header = header, double_precision = false,
                            ## quality control:
                            quality_control = quality_control, MAF = MAF, missing_value = missing_value,
                            ## others:
                            center = center, starting_value = starting_value)
        
        push!(layers[1].data, geno)

        #set the name of genotypes
        if is_fully_connected_12
            layers[1].data[1].name = layers[1].layer_name
            printstyled("In fully connected network, the name of genotypes (1st layer) is set to $(layers[1].layer_name).\n",bold=false,color=:green)
        else #partially connected 1->2 layer, set the name of genotypes for each omics
            layers[1].data[i].name = layers[1].layer_name * string(i) #e.g., "geno1", "geno2", "geno3"
            printstyled("In partially connected network, the name of $(i)th genotype (1st layer) is set to $(layers[1].data[i].name).\n",bold=false,color=:green)
        end
        push!(genotype_name_list, layers[1].data[i].name)
    end

    printstyled("End reading genotypes.......\n \n",bold=false,color=:green)
    printstyled("********************************************************** \n",bold=false,color=:green)

    ############################################################################
    # Step2. read omics in 2nd layer
    ############################################################################
    if layers[2].data != [] #not empty, user should re-run layers function to reset data stored in layer 2
        error("There is data already stored in layer 2, please re-run the code to reset layers.")
    end

    printstyled("Start reading omics .......\n",bold=false,color=:green)
    #arguments for data
    file_path       = layers[2].data_path
    separator       = layers[2].separator
    header          = layers[2].header
    quality_control = layers[2].quality_control = false #QC=false -> missing omics will not be replaced by its mean
    MAF             = layers[2].MAF             = false
    center          = layers[2].center          = false
    missing_value   = layers[2].missing_value
    
    #if file_path is a list of string (i.e., multiple omics files), error
    if typeof(file_path) == Vector{String}
        error("only one omics file is allowed in the 2nd layer.")
    end

    println()
    printstyled("Reading omics from $file_path \n",bold=false,color=:green)
    printstyled("Note that `center`, `quality_control`, and `MAF` in Layer() function are set to false here forcefully, so please use pre-processed omics data \n",bold=false,color=:green)

    #arguments for equations
    omics_name = equations[1].omics_name
    if omics_name==false
        error("please provide omics name in the 1st equation.")
    end

    #arguments for mme for 2->3 layer
    G = equations[2].G
    method = equations[2].method
    if method == "GBLUP"
        error("GBLUP is not supported from 2->3 layers. Try BayesC, RR-BLUP, etc..") #because when omics need imputation, GRM should be calculated in every MCMC iteration
    end
    printstyled("Note that `constraint` is set to false here forcefully \n",bold=false,color=:green)
    Pi                   = equations[2].Pi
    estimatePi           = equations[2].estimatePi
    G_is_marker_variance = equations[2].G_is_marker_variance #G is omics variance
    df_G                   = equations[2].df_G
    estimate_variance_G    = equations[2].estimate_variance_G
    estimate_scale_G       = equations[2].estimate_scale_G
    constraint_G           = equations[2].constraint_G

    if equations[2].starting_value != false
        error("user-defined starting value is not supported from 2->3 layer.")
    end

    # omics is of type "genotypes"
    omics = nnmm_get_omics(file_path, G;
                    ## omics name:
                    omics_name = omics_name,
                    ## method:
                    method = method, Pi = Pi, estimatePi = estimatePi, 
                    ## variance:
                    G_is_marker_variance = G_is_marker_variance, df = df_G,
                    estimate_variance = estimate_variance_G, estimate_scale = estimate_scale_G,
                    constraint = constraint_G, 
                    ## format:
                    separator = separator, header = true,
                    ## quality control:
                    missing_value = missing_value)
    #reset some parameters for omics
    # omics.alleleFreq = false
    # omics.sum2pq     = false


    # #check if data in 2nd layer are all missing
    # if sum(ismissing.(layers[2].data.genotypes)) == size(layers[2].data.genotypes,1) * size(layers[2].data.genotypes,2)
    #     printstyled("  - no omics data!! 2nd layer are all hidden nodes.\n",bold=false,color=:green)
    # end


    #save omics data in 2nd layer
    push!(layers[2].data, omics)
    
    printstyled("End reading omics.......\n \n",bold=false,color=:green)

    #for partial connected nn, #genotype files must = #omics in the 2nd layer
    if equations[1].partial_connect_structure == true
        if length(layers[1].data) != length(omics_name)
            error("For partial connected 1->2 layer, the number of genotype files must be equal to the number of omics in the 2nd layer.
                   number of genotype files: $(length(layers[1].data)),
                   number of omics in the 2nd layer: $(length(omics_name))")
        end
    end

    ############################################################################
    #Step2.5. read phenotypes
    ############################################################################
    if layers[3].data != [] #not empty, user should re-run layers function to reset data stored in layer 3
        error("There is data already stored in layer 3, please re-run the code to reset layers.")
    end

    printstyled("********************************************************** \n",bold=false,color=:green)
    printstyled("Start reading phenotypes.......\n",bold=false,color=:green)

    file_path       = layers[3].data_path
    separator       = layers[3].separator
    header          = layers[3].header
    missing_value   = layers[3].missing_value

    println()
    printstyled("Reading phenotypes from $file_path \n",bold=false,color=:green)

    pheno = read_phenotypes(file_path;
                            separator=separator,
                            header=header,
                            missing_value=missing_value)
    phenoID = equations[2].phenotype_name
    println("phenotype name: ",phenoID)
    if length(phenoID) > 1
        error("Multiple phenotypes in 3rd layer are not supported for now.")
    end
    
    layers[3].data=pheno #df with ID and other covariates/random effects data
    # #multi-trait is not supported for now
    # if layers[3].data.nPheno > 1
    #     error("Multi-trait in 3rd layer is not supported for now.")
    # end
    
    #show first 5 rows of phenotypes
    if size(layers[3].data,1)>5
        printstyled("First 5 rows of phenotypes (and other covariates/random effects data): \n",bold=false,color=:green)
        display(layers[3].data[1:5,:])
    end

    printstyled("End reading phenotypes.......\n \n",bold=false,color=:green)


    ############################################################################
    #Step3. check input equations
    ############################################################################
    printstyled("********************************************************** \n",bold=false,color=:green)
    printstyled("Start checking input data....... \n",bold=false,color=:green)
    println()
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

    # #check partial_connect_structure
    # is_fully_connected_12 = true #_12 means from 1st layer to 2nd layer
    # if equations[1].partial_connect_structure != false
    #     #make a marker effect model is used for partial connected 1->2 layer
    #     if equations[1].method=="GBLUP"
    #         error("The partial_connect_structure is not allowed for GBLUP model.")
    #     end
    #     #check if all elements are either 0 or 1, change to sparse matrix
    #     if all(equations[1].partial_connect_structure .== 0 .|| equations[1].partial_connect_structure .== 1)
    #         #change to sparse boolean matrix
    #         equations[1].partial_connect_structure = sparse(equations[1].partial_connect_structure .== 1) #change to sparse Boolean matrix
    #     else
    #         error("Elements in partial_connect_structure should be either 0 or 1.")
    #     end
    #     #check the dimention
    #     # - number of snps in the 1st layer = layers[1].data.nMarkers
    #     # - number of omics in the 2nd layer = layers[2].data.nMarkers
    #     if (layers[1].data.nMarkers,layers[2].data.nMarkers) != size(equations[1].partial_connect_structure)
    #         error("""
    #         Dimension mismatch in `partial_connect_structure` for layer 1 → 2.
            
    #         Current dimension: $(size(equations[1].partial_connect_structure))
    #         Expected dimension: (#nodes_layer1 × #nodes_layer2)
            
    #         Where:
    #           • #nodes_layer1 = $(layers[1].data.nMarkers)
    #           • #nodes_layer2 = $(layers[2].data.nMarkers)
            
    #         Please adjust `partial_connect_structure` to match the expected size.
    #         """)
            
    #     end
    #     is_fully_connected_12 = false
    # end

    # is_fully_connected_23 = true #_23 means from 2nd layer to 3rd layer, user-defined partial network is not allowed

    if equations[2].partial_connect_structure != false
        error("The partial_connect_structure is currently not allowed from 2->3 layers.")
    end

    #check nn activation function / user-defined function
    if isa(equations[2].activation_function, Function)  #e.g., user-defined function=Pig_Growth_Model(o1, o2, o3)
        num_args = first(methods(equations[2].activation_function)).nargs-1 #number of arguments in activation_function, should = number of nodes in the 2nd layer
        if num_args != layers[2].data.nMarkers #number of omics in the 2nd layer = layers[2].data.nMarkers
            error("The number of arguments in user-defined activation_function ≠ number of nodes in the 2nd layer.")
        end
        is_nn_activation_fcn = false #is user-defined function
    elseif equations[2].activation_function in ["tanh","sigmoid","relu","leakyrelu","linear"]
        is_nn_activation_fcn = true #is nn activation function
    else
        error("invalid activation_function.")
    end

    printstyled("End checking input data....... \n \n",bold=false,color=:green)

    ############################################################################
    #Step4. print NNMM info
    ############################################################################
    printstyled("******************** NNMM information ******************** \n",bold=false,color=:green)
    if is_fully_connected_12 #fully connected 1->2 layer
        printstyled(" - From 1st->2nd layers:             fully connected network. \n",bold=false,color=:green)
        printstyled("   - Number of nodes:     $(layers[1].data[1].nMarkers). \n",bold=false,color=:green)
    elseif !is_fully_connected_12 #partially connected 1->2 layer
        printstyled(" - From 1st->2nd layers:             partially connected network. \n",bold=false,color=:green)
        for i in 1:length(layers[1].data)
            printstyled("   - Number of nodes in $(i)th partial network:     $(layers[1].data[i].nMarkers). \n",bold=false,color=:green)
        end
    else
        error("error")
    end
    
    printstyled(" - From 2nd->3rd layers: \n",bold=false,color=:green)
    printstyled("   - Number of nodes:     $(layers[2].data[1].nFeatures). \n",bold=false,color=:green)

    #print nn activation function / user-defined function
    if is_nn_activation_fcn  #nn activation function such as tanh, sigmoid
        printstyled(" - Activation function in 2nd layer: $(equations[2].activation_function).\n",bold=false,color=:green)
        printstyled(" - Sample missing data in 2nd layer: Hamiltonian Monte Carlo. \n",bold=false,color=:green)
    elseif !is_nn_activation_fcn #user-defined function. e.g, Pig_Growth_Model(o1, o2, o3)
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
    if is_fully_connected_12
        equation_12 = ""
        _, equ_right = strip.(split(equations[1].equation,"=")) # "intercept + x1 + geno" 
        for i in layers[2].data[1].featureID #the name of omics
            equation_12 = equation_12 * i * "=" * equ_right * ";"
        end
    else #partial connected 1->2 layer
        equation_12 = ""
        for i in 1:length(layers[1].data)
            _, equ_right = strip.(split(equations[1].equation,"=")) # "intercept + x1 + geno"
            #in equ_right, replace layer name (e.g., "geno") with the name of the i-th genotype (e.g., "geno1")
            equ_right = replace(equ_right, layers[1].layer_name => layers[1].data[i].name)
            omics_name = layers[2].data[1].featureID[i]
            equation_12 = equation_12 * omics_name * "=" * equ_right * ";"
        end
    end
    equation_12 = equation_12[1:(end-1)] #remove the ";" at the end
    equations[1].equation = equation_12 #update the equation
    
    #print equations
    printstyled(" - Model equation 1->2:\n", bold=false, color=:green)
    for eq in split(equation_12, ';')
        printstyled("    $eq\n", bold=false, color=:green)
    end

    ############################################################################
    # Step6. build model equation (2->3): replace "phenotypes" with trait names
    #  - old: phenotypes=intercept+omics
    #  - new: y1=intercept+omics
    ############################################################################
    equation_23 = ""
    _, equ_right = strip.(split(equations[2].equation,"=")) # "intercept + omics"
   
    for i in phenoID
        equation_23 = equation_23 * i * "=" * equ_right * ";"
    end
    equation_23 = equation_23[1:(end-1)] #remove the ";" at the end
    equations[2].equation = equation_23 #update the equation
    
    #print equations
    printstyled(" - Model equation 2->3:\n", bold=false, color=:green)
    for eq in split(equation_23, ';')
        printstyled("    $eq\n", bold=false, color=:green)
    end

    mme_all = []
    ############################################################################
    # 1->2: All model terms (will be added to MME)
    # same as build_model()
    ############################################################################
    printstyled("********************************************************** \n",bold=false,color=:green)
    printstyled("Start build mme for 1->2....... \n",bold=false,color=:green)

    if equations[1].R != false && !isposdef(map(AbstractFloat,equations[1].R))
      error("The covariance matrix is not positive definite.")
    end
    if !(typeof(equations[1].equation)<:AbstractString) || equations[1].equation==""
      error("Model equations are wrong.\n
      To find an example, type ?build_model and press enter.\n")
    end
    if equations[1].estimate_scale_R != false
      error("estimate scale for residual variance is not supported now.")
    end

    model_equations=equations[1].equation
    #e.g., ""y2 = A+B+A*B""
    modelVec   = [strip(i) for i in split(model_equations,[';','\n'],keepempty=false)]
    nModels    = size(modelVec,1)
    if equations[1].R != false && size(equations[1].R,1) != nModels
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
    genotypes = [] #note that the genotypei is reference from layers[1].data
    whichterm = 1
    for term in modelTerms
        term_symbol = Symbol(split(term.trmStr,":")[end])
        if term_symbol ∈ Symbol.(genotype_name_list)
          term.random_type = "genotypes"
          if is_fully_connected_12
            genotypei = layers[1].data[1]
          else
            genotypei = layers[1].data[findfirst(x->x==term_symbol, Symbol.(genotype_name_list))]
          end
          if genotypei.name ∉ map(x->x.name, genotypes) #only save unique genotype
            genotypei.ntraits = is_fully_connected_12 ? nModels : 1
            genotypei.trait_names = is_fully_connected_12 ? string.(lhsVec) : term.iTrait
            # mega_trait
            # if nModels != 1 
            #   genotypei.G.df = genotypei.G.df + nModels
            # end
            if is_fully_connected_12 && (genotypei.G.val != false || genotypei.genetic_variance.val != false)
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
    R=equations[1].R
    df_R=equations[1].df_R
    estimate_variance_R=equations[1].estimate_variance_R
    estimate_scale_R=equations[1].estimate_scale_R
    constraint_R=equations[1].constraint_R
    if constraint_R == false
        error("constraint_R must be true for mega_trait.")
    end
    # if nModels == 1
        scale_R = R*(df_R - 2)/df_R
        df_R    = df_R 
    # else
    #     scale_R = R*(df_R - 1)
    #     df_R    = df_R + nModels
    # end

    #initialize mme
    mme1 = MME(nModels,modelVec,modelTerms,dict,lhsVec, 
    Variance(R==false ? R : Float32.(R), #val
            Float32(df_R),              #df
            R==false ? R : scale_R,     #scale
            estimate_variance_R, estimate_scale_R, constraint_R))

    mme1.M = genotypes #add genotypes into mme
    
    mme1.is_fully_connected   = is_fully_connected_12
    mme1.is_activation_fcn    = is_nn_activation_fcn
    mme1.nonlinear_function   = isa(equations[2].activation_function, Function) ? equations[2].activation_function : nnbayes_activation(equations[2].activation_function)
    # mme.yobs_name=mme.lhsVec
    mme1.latent_traits = true

    ############################################################################
    # 1->2: set random
    ############################################################################
    if equations[1].random != false
        for r in equations[1].random
            @assert haskey(r, :name) "Each spec must have a :name"
            name = r.name

            if haskey(r, :pedigree)
                # pedigree path
                G = get(r, :G, false)
                set_random(mme1, name, r.pedigree, G)
                printstyled("set random term for $(name) with pedigree \n",bold=false,color=:green)
            elseif haskey(r, :Vinv) 
                Vinv = Matrix{Float64}(r.Vinv)
                G = get(r, :G, false)
                index = get(r, :index, String[])
                set_random(mme1, name, G; Vinv=Vinv, names=index)
                printstyled("set random term for $(name) with user-provided relationship matrix (inversed) \n",bold=false,color=:green) 
            else
                # simple path
                set_random(mme1, name)
                printstyled("set random term for $(name) \n",bold=false,color=:green)
            end
    end

    ############################################################################
    # 1->2: set covariate
    ############################################################################
    if equations[1].covariate !=false
        for covariate_term in equations[1].covariate
            if !(covariate_term isa String)
                error("double-check covariate term input")
            end
            set_covariate(mme1,covariate_term)
            printstyled("set covariate term for $(covariate_term) \n",bold=false,color=:green)
        end
    end

    push!(mme_all,mme1)
    printstyled("End build mme for layer 1->2....... \n \n",bold=false,color=:green)
    
    ############################################################################
    # 2->3: All model terms (will be added to MME)
    # same as build_model()
    ############################################################################
    printstyled("********************************************************** \n",bold=false,color=:green)
    printstyled("Start build mme for layer 2->3....... \n",bold=false,color=:green)
    
    is_fully_connected_23 = true #must be true now

    if equations[2].R != false && !isposdef(map(AbstractFloat,equations[2].R))
      error("The covariance matrix is not positive definite.")
    end
    if !(typeof(equations[2].equation)<:AbstractString) || equations[2].equation==""
      error("Model equations are wrong.\n
      To find an example, type ?build_model and press enter.\n")
    end
    if equations[2].estimate_scale_R != false
      error("estimate scale for residual variance is not supported now.")
    end
    
    model_equations=equations[2].equation
    #e.g., ""y2 = A+B+A*B""
    modelVec   = [strip(i) for i in split(model_equations,[';','\n'],keepempty=false)]
    nModels    = size(modelVec,1)

    if equations[2].R != false && size(equations[2].R,1) != nModels
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
    omics = []
    whichterm = 1
    for term in modelTerms
        term_symbol = Symbol(split(term.trmStr,":")[end])
        if term_symbol == Symbol(layers[2].layer_name) #e.g., :omics
          term.random_type = "genotypes" #omics
          omicsi = layers[2].data[1]
          omicsi.name = string(term_symbol)
          if omicsi.name ∉ map(x->x.name, omics) #only save unique omics
            omicsi.ntraits = nModels
            omicsi.trait_names = string.(lhsVec)
            if nModels != 1
                omicsi.G.df = omicsi.G.df + nModels
            end
            if is_fully_connected_23 && (omicsi.G.val != false || omicsi.genetic_variance.val != false)
              if size(omicsi.G.val,1) != nModels && size(omicsi.genetic_variance.val,1) != nModels
                error("The genomic covariance matrix is not a ",nModels," by ",nModels," matrix.")
              end
            end
            push!(omics,omicsi)
          end
        end
    end

    #create mme with omics
    filter!(x->x.random_type != "genotypes",modelTerms) #remove "genotypes" from modelTerms
    filter!(x->x[2].random_type != "genotypes",dict)    #remove "genotypes" from dict
    
    #set scale and df for residual variance
    R=equations[2].R
    df_R=equations[2].df_R
    estimate_variance_R=equations[2].estimate_variance_R
    estimate_scale_R=equations[2].estimate_scale_R
    constraint_R=equations[2].constraint_R
    if nModels == 1
        scale_R = R*(df_R - 2)/df_R
        df_R    = df_R 
    else
        scale_R = R*(df_R - 1)
        df_R    = df_R + nModels
    end

    #initialize mme
    mme2 = MME(nModels,modelVec,modelTerms,dict,lhsVec, 
    Variance(R==false ? R : Float32.(R), #val
            Float32(df_R),              #df
            R==false ? R : scale_R,     #scale
            estimate_variance_R, estimate_scale_R, constraint_R))

    mme2.M = omics #add omics into mme
    mme2.is_fully_connected   = is_fully_connected_23 # =true
    

    mme2.is_activation_fcn    = false #this is only for 1->2 layer now
    mme2.nonlinear_function   = false #this is only for 1->2 layer now
    # mme.yobs_name=mme.lhsVec
    mme2.latent_traits = false

    ############################################################################
    # 2->3: set random in all models 
    ############################################################################
    if equations[2].random != false
        for r in equations[2].random
            @assert haskey(r, :name) "Each spec must have a :name"
            name = r.name

            if haskey(r, :pedigree)
                # pedigree path
                G = get(r, :G, false)
                set_random(mme2, name, r.pedigree, G)
                printstyled("set random term for $(name) with pedigree \n",bold=false,color=:green)
            elseif haskey(r, :Vinv) 
                Vinv = Matrix{Float64}(r.Vinv)
                G = get(r, :G, false)
                index = get(r, :index, String[])
                set_random(mme2, name, G; Vinv=Vinv, names=index)
                printstyled("set random term for $(name) with user-provided relationship matrix (inversed) \n",bold=false,color=:green) 
            else
                # simple path
                set_random(mme2, name)
                printstyled("set random term for $(name) \n",bold=false,color=:green)
            end
    end

    ############################################################################
    # 2->3: set covariate in all models 
    ############################################################################
    if equations[2].covariate !=false
        for covariate_term in equations[2].covariate
            if !(covariate_term isa String)
                error("double-check covariate term input")
            end
            set_covariate(mme2,covariate_term)
            printstyled("set covariate term for $(covariate_term) \n",bold=false,color=:green)
        end
    end

    push!(mme_all,mme2)
    
    printstyled("End build mme for 2->3....... \n \n",bold=false,color=:green)
    

    ############################################################################
    # run MCMC (from runMCMC())
    ############################################################################
    printstyled("********************************************************** \n",bold=false,color=:green)
    printstyled("Start prepare MCMC....... \n",bold=false,color=:green)
    # 1->2
    mme_all[1].yobs_name=equations[2].phenotype_name[1] #yobs is phenotype in 3rd layer

    
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
    heterogeneous_residuals=false
    single_step_analysis=false
    fitting_J_vector=false
    missing_phenotypes=false #NN-MM -> missing hidden nodes will be sampled
    RRM=false
    fast_blocks=false

    #1->2
    mme_all[1].MCMCinfo = MCMCinfo(heterogeneous_residuals,
                   chain_length,burnin,output_samples_frequency,
                   printout_model_info,printout_frequency, single_step_analysis,
                   fitting_J_vector,missing_phenotypes,
                   update_priors_frequency,outputEBV,output_heritability,prediction_equation,
                   seed,double_precision,output_folder,RRM,fast_blocks)

    #2->3
    missing_phenotypes=true # Bayesian Alphabet
    mme_all[2].MCMCinfo = MCMCinfo(heterogeneous_residuals,
                            chain_length,burnin,output_samples_frequency,
                            printout_model_info,printout_frequency, single_step_analysis,
                            fitting_J_vector,
                            missing_phenotypes, ##NN-MM -> missing hidden nodes will be sampled
                            update_priors_frequency,outputEBV,output_heritability,prediction_equation,
                            seed,double_precision,output_folder,RRM,fast_blocks)

    ############################################################################
    # Check 1)Arguments; 2)Input Pedigree,Genotype,Phenotypes,
    #       3)output individual IDs; 4)Priors  5)Prediction equation
    #       in the Whole Dataset.
    ############################################################################
    #1->2
    printstyled("check errors args for 1->2 \n",bold=false,color=:green)
    errors_args(mme_all[1])       #check errors in function arguments
    
    #2->3
    printstyled("check errors args for 2->3 \n",bold=false,color=:green)
    errors_args(mme_all[2])       #check errors in function arguments

    printstyled("check pedigree genotypes phenotypes \n",bold=false,color=:green)
    #merge phenotypes with omics data
    layers[2].data[1].data=leftjoin(layers[2].data[1].data, layers[3].data, on=:ID)

    # Find phenotype rows that have no match in omics
    unmatched_df2 = antijoin(layers[3].data, layers[2].data[1].data, on=:ID)  # rows in phenotype with no key in omics
    if nrow(unmatched_df2) > 0
        @warn "IDs in phenotypes data not matched by omics data" unmatched_df2
    end
    #1->2
    layers[2].data[1].data=check_pedigree_genotypes_phenotypes(mme_all[1],layers[2].data[1].data,false) #pedigree = false for none single-step
    #2->3
    layers[2].data[1].data=check_pedigree_genotypes_phenotypes(mme_all[2],layers[2].data[1].data,false) #pedigree = false for none single-step

    #initialize missing omics data
    printstyled("initialize missing omics data \n",bold=false,color=:green)
    #initiliza missing omics data  (after check_pedigree_genotypes_phenotypes() because non-genotyped inds are removed)
    nnlmm_initialize_missing_with_mean(mme_all[1],layers[2].data[1].data)

    #1->2
    printstyled("1->2 prediction setup: \n",bold=false,color=:green)
    prediction_setup(mme_all[1])  #set prediction equation, defaulting to genetic values
    check_outputID(mme_all[1])    #check individual of interest for prediction
    println()

    #2->3
    printstyled("2->3 prediction setup: \n",bold=false,color=:green)
    prediction_setup(mme_all[2])  #set prediction equation, defaulting to genetic values
    check_outputID(mme_all[2])    #check individual of interest for prediction
    println()


    #1->2
    printstyled("1->2 make_dataframes: \n",bold=false,color=:green)
    df_whole1,train_index1 = make_dataframes(layers[2].data[1].data,mme_all[1])
    println()
    
    #2->3
    printstyled("2->3 make_dataframes: \n",bold=false,color=:green)
    df_whole2,train_index2 = make_dataframes(layers[2].data[1].data,mme_all[2])
    println()


    #1->2
    printstyled("1->2 set_default_priors_for_variance_components: \n",bold=false,color=:green)
    set_default_priors_for_variance_components(mme_all[1],df_whole1)  #check priors (set default priors)
    set_marker_hyperparameters_variances_and_pi(mme_all[1])
    println()
    
    #2->3
    printstyled("2->3 set_default_priors_for_variance_components: \n",bold=false,color=:green)
    set_default_priors_for_variance_components(mme_all[2],df_whole2)  #check priors (set default priors)
    set_marker_hyperparameters_variances_and_pi(mme_all[2])
    println()

    ############################################################################
    # Adhoc functions
    ############################################################################
    # Double Precision
    if double_precision == true
        for mme in mme_all
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
    end


    # NNBayes: modify parameters for partial connected NN
    if equations[1].partial_connect_structure==true
        nnbayes_partial_para_modify2(mme_all[1])
    end
    ############################################################################
    #Make incidence matrices and genotype covariates for training observations
    #and individuals of interest
    ############################################################################
    #make incidence matrices (non-genomic effects) (after SSBRrun for ϵ & J)
    layers[2].data[1].data=make_incidence_matrices(mme_all[1],df_whole1,train_index1) #`layers[2].data[1].data` is `df`
    df2=make_incidence_matrices(mme_all[2],df_whole2,train_index2)

    #align genotypes with 1) phenotypes IDs; 2) output IDs.
    align_genotypes(mme_all[1],output_heritability,single_step_analysis)
    # initiate Mixed Model Equations and check starting values
    #1->2
    init_mixed_model_equations(mme_all[1], layers[2].data[1].data, false) #(mme,df,starting_value)

    #2->3
    #align omics with output IDs.
    align_omics(mme_all[2],output_heritability,single_step_analysis)
    #align genotypes with phenotypes IDs;
    align_transformed_omics_with_phenotypes(mme_all[2],mme_all[1].nonlinear_function)
    # initiate Mixed Model Equations and check starting values
    init_mixed_model_equations(mme_all[2], df2, false) #(mme,df,starting_value)


    ############################################################################
    # MCMC
    ############################################################################
    #1->2
    printstyled("Describe layer 1->2: \n",bold=false,color=:green)
    describe(mme_all[1])

    #2->3
    printstyled("Describe layer 2->3: \n",bold=false,color=:green)
    describe(mme_all[2])

    t=nnmm_MCMC_BayesianAlphabet(mme_all[1],layers[2].data[1].data, mme_all[2], df2) #mme1,df1,mme2,df2
    return t

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
end
end