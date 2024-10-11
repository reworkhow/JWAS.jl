function nnmm_runMCMC(mme::MME,df;
                #Data
                heterogeneous_residuals           = false,
                #MCMC
                chain_length::Integer             = 100,
                starting_value                    = false,
                burnin::Integer                   = 0,
                output_samples_frequency::Integer = chain_length>1000 ? div(chain_length,1000) : 1,
                update_priors_frequency::Integer  = 0,
                #Methods
                single_step_analysis            = false, #parameters for single-step analysis
                pedigree                        = false, #parameters for single-step analysis
                fitting_J_vector                = true,  #parameters for single-step analysis
                causal_structure                = false,
                missing_phenotypes              = false,
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
    mme.output=MCMC_BayesianAlphabet(mme,df)

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