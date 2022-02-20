################################################################################
# Check Input Data before running MCMC for
# 1) arguments
# 2) output IDs
# 3) genotypes/pedigree/phenotypes
# 4) priors
# 6) starting values
################################################################################
function errors_args(mme)
    if mme.mmePos != 1
      error("Please build your model again using the function build_model().")
    end

    if mme.MCMCinfo.output_samples_frequency <= 0
        error("output_samples_frequency should be an integer > 0.")
    end

    if mme.M != 0
        for Mi in mme.M
            if !(Mi.method in ["BayesL","BayesC","BayesB","BayesA","RR-BLUP","GBLUP"])
                error(Mi.method," is not available in JWAS. Please read the documentation.")
            end

            if Mi.method in ["RR-BLUP","BayesL","GBLUP","BayesA"]
                if Mi.π != false
                    printstyled(Mi.method," runs with π = false.\n",bold=false,color=:red)
                    Mi.π = false
                elseif Mi.estimatePi == true
                    printstyled(Mi.method," runs with estimatePi = false.\n",bold=false,color=:red)
                    Mi.estimatePi = false
                end
            end
            if Mi.method == "BayesA"
                Mi.method = "BayesB"
                println("BayesA is equivalent to BayesB with known π=0. BayesB with known π=0 runs.")
            end
            if Mi.method == "GBLUP"
                if Mi.genetic_variance == false && Mi.G != false
                    error("Please provide values for the genetic variance for GBLUP analysis")
                end
                if mme.MCMCinfo.single_step_analysis == true
                    error("SSGBLUP is not available")
                end
            end
            if mme.nModels > 1 && Mi.π != 0.0
                if round(sum(values(Mi.π)),digits=2)!=1.0
                  error("Summation of probabilities of Pi is not equal to one.")
                end
                if typeof(Mi.π) <: Number
                    error("Pi cannot be a number in multi-trait analysis.")
                end
            end
        end
    end
    if mme.MCMCinfo.single_step_analysis == true && mme.M == 0
        error("Genomic information is required for single-step analysis.")
    end
    if mme.causal_structure != false && mme.nModels == 1
        error("Causal strutures are only allowed in multi-trait analysis")
    end
end

function check_outputID(mme)
    #***************************************************************************
    #Set default output IDs
    #
    #Genotyped individuals are usaully not many, and are used in GWAS (complete
    #and incomplete), thus are used as default output_ID if not provided
    #***************************************************************************
    if mme.MCMCinfo.outputEBV == true
        if mme.output_ID == false
            if mme.M != 0   #genomic analysis: all genotyped IDs
                mme.output_ID = copy(mme.M[1].obsID)
            else            #non-genomic analyis
                if mme.pedTrmVec != 0 #all pedigreed IDs in Pedigree-based BLUP
                    mme.output_ID = copy(mme.ped.IDs)
                else                  #no EBV in non-genetic analysis
                    mme.MCMCinfo.outputEBV = false
                end
            end
        end
    else
        mme.output_ID = false
    end
    #***************************************************************************
    #Set ouput IDs in heritability (h²) estimation
    #***************************************************************************
    if mme.MCMCinfo.output_heritability == true && mme.M != 0
        mme.MCMCinfo.outputEBV == true
        if mme.MCMCinfo.single_step_analysis == false
            mme.output_ID = copy(mme.M[1].obsID)
        else
            mme.output_ID = copy(mme.ped.IDs)
        end
    end
    #***************************************************************************
    #Remove output IDs not found in genotypes or pedigree
    #***************************************************************************
    single_step_analysis = mme.MCMCinfo.single_step_analysis
    if single_step_analysis == false && mme.M != 0 #complete genomic data
        if mme.output_ID!=0 && !issubset(mme.output_ID,mme.M[1].obsID)
            printstyled("Testing individuals are not a subset of ",
            "genotyped individuals (complete genomic data,non-single-step). ",
            "Only output EBV for tesing individuals with genotypes.\n",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,mme.M[1].obsID)
        end
    elseif mme.ped != false #1)incomplete genomic data 2)PBLUP
        pedID = map(string,collect(keys(mme.ped.idMap)))
        if mme.output_ID!=0 && !issubset(mme.output_ID,pedID)
            printstyled("Testing individuals are not a subset of ",
            "individuals in pedigree (incomplete genomic data (single-step) or PBLUP). ",
            "Only output EBV for tesing individuals in the pedigree.",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,pedID)
        end
    end
end

function check_pedigree_genotypes_phenotypes(mme,df,pedigree)
    #1)Pedigree
    if mme.ped != false || pedigree != false
        printstyled("Checking pedigree...\n" ,bold=false,color=:green)
        if mme.ped == false #check whether mme.ped exists
            mme.ped = deepcopy(pedigree)
        end
    end
    #2)Genotypes
    if mme.M != false
        printstyled("Checking genotypes...\n" ,bold=false,color=:green)
        for Mi in mme.M
            if Mi.obsID != mme.M[1].obsID
                error("genotypic information is not provided for same individuals")
            end
            if mme.ped != 0
                pedID=map(string,collect(keys(mme.ped.idMap)))
                if !issubset(Mi.obsID,pedID)
                    error("Not all genotyped individuals are found in pedigree!")
                end
            end
        end
        if mme.MCMCinfo.single_step_analysis == true && length(mme.M) != 1
            error("Now only one genomic category is allowed in single-step analysis.")
        end
    end
    #3)Phenotypes
    printstyled("Checking phenotypes...\n" ,bold=false,color=:green)
    printstyled("Individual IDs (strings) are provided in the first column of the phenotypic data.\n" ,bold=false,color=:green)
    ##remove records according to genotype and pedigree information
    phenoID = df[!,1] = map(string,strip.(map(string,df[!,1]))) #make IDs stripped string
    if mme.M != false && mme.MCMCinfo.single_step_analysis == false #complete genomic data
        if !issubset(phenoID,mme.M[1].obsID)
            index = [phenoID[i] in mme.M[1].obsID for i=1:length(phenoID)]
            df    = df[index,:]
            printstyled("In this complete genomic data (non-single-step) analyis, ",
            length(index)-sum(index)," phenotyped individuals are not genotyped. ",
            "These are removed from the analysis.\n",bold=false,color=:red)
        end
    end
    if mme.ped != false #1)incomplete genomic data 2)PBLUP or 3)complete genomic data with polygenic effect
        pedID = map(string,collect(keys(mme.ped.idMap)))
        if !issubset(phenoID,pedID)
            index = [phenoID[i] in pedID for i=1:length(phenoID)]
            df    = df[index,:]
            if mme.MCMCinfo.censored_trait != false
                mme.MCMCinfo.censored_trait = mme.MCMCinfo.censored_trait[index]
            end
            printstyled("In this incomplete genomic data (single-step) or PBLUP analysis, ",
            length(index)-sum(index)," phenotyped individuals are not included in the pedigree. ",
            "These are removed from the analysis.\n",bold=false,color=:red)
        end
        if mme.M != false
            if issubset(phenoID,mme.M[1].obsID) && mme.MCMCinfo.single_step_analysis == true
                error("All phenotyped individuals arer genotyped. Do not run single step analysis.")
            end
        end
    end
    rename!(df,strip.(names(df)))
    return df
end

function set_default_priors_for_variance_components(mme,df)
  myvar     = [var(filter(isfinite,skipmissing(df[!,mme.lhsVec[i]]))) for i=1:size(mme.lhsVec,1)]
  phenovar  = diagm(0=>myvar)
  h2        = 0.5

  genetic_random_count    = (mme.M!=0 ? length(mme.M) : 0)
  nongenetic_random_count = 1
  for random_term in mme.rndTrmVec
    if random_term.randomType == "A"
        genetic_random_count += 1
    elseif split(random_term.term_array[1],":")[2] != "ϵ"
        nongenetic_random_count +=1
    end
  end
  if genetic_random_count != 0
      varg = phenovar*h2/genetic_random_count
  end
  vare = phenovar*h2/nongenetic_random_count

  #genetic variance or marker effect variance
  if mme.M!=0
    for Mi in mme.M
      if Mi.G == false && Mi.genetic_variance == false
          printstyled("Prior information for genomic variance is not provided and is generated from the data.\n",bold=false,color=:green)
          if mme.nModels==1
              Mi.genetic_variance = varg[1,1]
          elseif mme.nModels>1
              Mi.genetic_variance = varg
          end #mme.M.G and its scale parameter will be reset in function set_marker_hyperparameters_variances_and_pi
      end
    end
  end
  #residual effects
  if mme.nModels==1 && isposdef(mme.R) == false #single-trait
    printstyled("Prior information for residual variance is not provided and is generated from the data.\n",bold=false,color=:green)
    is_categorical_trait = mme.traits_type==["categorical"]
    mme.R = mme.ROld = is_categorical_trait ? 1.0 : vare[1,1]  #residual variance known to be 1.0 in single trait categorical analysis
    mme.scaleR = mme.R*(mme.df.residual-2)/mme.df.residual
  elseif mme.nModels>1 && isposdef(mme.R) == false #multi-trait
    printstyled("Prior information for residual variance is not provided and is generated from the data.\n",bold=false,color=:green)
    mme.R = mme.ROld = vare
    mme.scaleR = mme.R*(mme.df.residual - mme.nModels - 1)
  end
  #random effects
  if length(mme.rndTrmVec) != 0
    for randomEffect in mme.rndTrmVec
      if isposdef(randomEffect.Gi) == false
        printstyled("Prior information for random effect variance is not provided and is generated from the data.\n",bold=false,color=:green)
        myvarout  = [split(i,":")[1] for i in randomEffect.term_array]
        myvarin   = string.(mme.lhsVec)
        Zdesign   = mkmat_incidence_factor(myvarout,myvarin)
        if randomEffect.randomType == "A"
            G = diagm(Zdesign*diag(varg))
        else
            G = diagm(Zdesign*diag(vare))
        end
        randomEffect.Gi = randomEffect.GiOld = randomEffect.GiNew = Symmetric(inv(G))
        randomEffect.scale = G*(randomEffect.df-length(randomEffect.term_array)-1)
        if randomEffect.randomType == "A"
          mme.Gi = randomEffect.Gi
          mme.scalePed = randomEffect.scale
        end
      end
    end
  end
end

#return a dataframe of all observations;
#missing values are used for some individuals of interest for which data is not available.
function make_dataframes(df,mme)
    #***************************************************************************
    #NN-Bayes Omics: individuals with all omics data but no yobs should be kept
    #                since the omics data can help better estimate marker effect
    #***************************************************************************
    if mme.nonlinear_function != false && mme.latent_traits != false
        lhsVec = [mme.yobs_name ; mme.lhsVec]  # [:y, :gene1, :gene2]
    else
        lhsVec = mme.lhsVec #reference, not copy
    end
    #***************************************************************************
    #Whole Data (training + individuals of interest)
    #***************************************************************************
    #expand phenotype dataframe to include individuals of interest (to make incidencee matrices)
    if !issubset(mme.output_ID,df[!,1]) && mme.output_ID != false
        #IDs for some individuals of interest are not included in the phenotypic data.
        #These are added as additional rows with missing values in df_whole.
        df_output = DataFrame(ID=setdiff(mme.output_ID,df[!,1]))
        rename!(df_output,"ID"=>names(df)[1])
        df_whole = outerjoin(df,df_output,on=names(df)[1])
    else
        df_whole = df
    end
    #***************************************************************************
    #Training data
    #(non-missing observations, only these ar used in mixed model equations)
    #***************************************************************************
    #remove individuals whose phenotypes are missing for all traits fitted in the model
    missingdf  = ismissing.(Matrix(df_whole[!,lhsVec]))
    allmissing = fill(true,length(lhsVec))
    train_index = Array{Int64,1}()
    for i in 1:size(missingdf,1)
        if missingdf[i,:] != allmissing
            push!(train_index,i)
        end
    end
    #***************************************************************************
    # Check missing levels in factors
    #***************************************************************************
    #check training data with model equations
    for term in mme.modelTerms
        for i = 1:term.nFactors
          if term.factors[i] != :intercept && any(ismissing,df_whole[train_index,term.factors[i]])
            printstyled("Missing values are found in column ",term.factors[i]," for some observations.",
            "Effects of this variable on such observations are considered as zeros. ",
            "It will be used in the estimation of this effect. Users may impute missing values before the analysis \n",
            bold=false,color=:red)
          end
        end
    end
    #check individual of interest (output IDs) with prediction equations
    for term in mme.MCMCinfo.prediction_equation
        term_value = mme.modelTermDict[term]
        for i = 1:term_value.nFactors
          if term_value.factors[i] != :intercept && any(ismissing,df_whole[!,term_value.factors[i]])
            printstyled("Missing values are found in column ",term_value.factors[i]," for some observations.",
            "Effects of this variable on such observations are considered as zeros. ",
            "It will be used in prediction. Users may impute missing values before the analysis. \n",
            bold=false,color=:red)
          end
        end
    end
    #***************************************************************************
    # Return
    #***************************************************************************
    printstyled("Phenotypes for ",length(train_index)," observations are used in the analysis.",
    "These individual IDs are saved in the file IDs_for_individuals_with_phenotypes.txt.\n",bold=false,color=:green)
    writedlm("IDs_for_individuals_with_phenotypes.txt",unique(df_whole[train_index,1]))
    mme.obsID = map(string,df_whole[train_index,1])
    return df_whole,train_index
end

"""
    make_incidence_matrices(mme,df_whole,train_index)

* (internal function) Make incidence matrices for effects involved in
calculation of EBV except marker covariates.
* Both incidence matrices for non-missing observations (used in mixed model equations)
and individuals of interest (output IDs) are obtained.
"""
function make_incidence_matrices(mme,df_whole,train_index)
    #***************************************************************************
    #Whole Data
    #***************************************************************************
    #Make incidence matrices X for each term (whole data)
    for term in mme.modelTerms
      getData(term,df_whole,mme)
      getX(term,mme)
    end
    #***************************************************************************
    #individuals of interest (output EBV)
    #***************************************************************************
    #check whehter all levels in output_X exist in training data
    #e.g, g1,g2,g6 in output but g6 doesn't exist in g1,g2
    for term in mme.MCMCinfo.prediction_equation
        Zout = mkmat_incidence_factor(string(mme.modelTermDict[term].iModel) .* mme.output_ID,
               vcat(Tuple([string(i) .* df_whole[!,1] for i=1:mme.nModels])...))
        mme.output_X[term] = Zout*mme.modelTermDict[term].X
    end
    #***************************************************************************
    #Training data
    #***************************************************************************
    for term in mme.modelTerms
        train_sel = [i in train_index for i=1:size(df_whole,1)]
        term.X = term.X[repeat(train_sel,mme.nModels),:]
    end
    df = df_whole[train_index,:]
    #require all levels for output id being observed in training data
    #for any fixed factor or i.i.d random factor
    for term in mme.modelTerms
        if term.nLevels > 1 && term.random_type in ["fixed","I"]
            train_effects = vec(sum(term.X,dims=1) .!= 0.0)
            if sum(train_effects) != length(train_effects) #where zero columns exist
                error("Some levels in $(term.trmStr) for individuals of interest are not found in training individuals (IDs with non-missing records).",
                      "You may delete those rows or replace those values with missing.")
            end
        end
    end
    return df
end

function init_mixed_model_equations(mme,df,starting_value)
    getMME(mme,df)
    ############################################################################
    #starting value for non-marker location parameters (sol) can be provided
    ############################################################################
    nsol = size(mme.mmeLhs,1)
    #nsol = size(mme.mmeLhs,1)+sum((Mi.method != "GBLUP" ? Mi.nMarkers : Mi.nObs) for Mi in mme.M)*mme.nModels
    if starting_value == false #no starting values
        mme.sol = zeros((mme.MCMCinfo.double_precision ? Float64 : Float32),nsol)
    else            #besure type is Float64
        printstyled("Starting values are provided. The order of starting values for location parameters\n",
        "should be the order of location parameters in the Mixed Model Equation for all traits (This can be\n",
        "obtained by getNames(model)).\n",bold=false,color=:green)
        if length(starting_value) != nsol
            error("length of starting values for non-marker location parameters is wrong.")
        end
        mme.sol = map((mme.MCMCinfo.double_precision ? Float64 : Float32),starting_value)
    end
    ############################################################################
    #starting value marker effects
    ############################################################################
    if mme.M != 0
        for Mi in mme.M
            nsol = Mi.method != "GBLUP" ? Mi.nMarkers : Mi.nObs
            if Mi.α == false
                Mi.α = zeros((mme.MCMCinfo.double_precision ? Float64 : Float32),nsol*Mi.ntraits)
            else
                if length(Mi.α) != nsol*Mi.ntraits
                    error("length of starting values for marker effects is wrong.")
                end
                Mi.α = map((mme.MCMCinfo.double_precision ? Float64 : Float32),Mi.α)
            end
            Mi.α = [Mi.α[(traiti-1)*nsol+1:traiti*nsol] for traiti = 1:Mi.ntraits]
        end
    end
end
