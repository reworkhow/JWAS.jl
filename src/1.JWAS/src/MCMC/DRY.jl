################################################################################
# Pre-Check
################################################################################
function errors_args(mme,methods)
    if mme.mmePos != 1
      error("Please build your model again using the function build_model().")
    end

    Pi         = mme.MCMCinfo.Pi
    estimatePi = mme.MCMCinfo.estimatePi
    if methods == "conventional (no markers)"
        if mme.M!=0
            error("Conventional analysis runs without genotypes!")
        elseif estimatePi == true
            error("conventional (no markers) analysis runs with estimatePi = false.")
        end
    elseif methods=="RR-BLUP"
        if mme.M == 0
            error("RR-BLUP runs with genotypes")
        elseif Pi != 0.0
            error("RR-BLUP runs with π=0.")
        elseif estimatePi == true
            error("RR-BLUP runs with estimatePi = false.")
        end
    elseif methods=="BayesC"
        if mme.M == 0
            error("BayesC runs with genotypes.")
        end
    elseif methods=="BayesB"
        if mme.M==0
            error("BayesB runs with genotypes.")
        end
    elseif methods=="BayesL"
        if mme.M == 0
            error("BayesL runs with genotypes.")
        elseif estimatePi == true
            error("BayesL runs with estimatePi = false.")
        end
    elseif methods=="GBLUP"
        if mme.M == 0
            error("GBLUP runs with genotypes.")
        elseif mme.M.genetic_variance == false
            error("Please provide values for the genetic variance for GBLUP analysis")
        elseif estimatePi == true
            error("GBLUP runs with estimatePi = false.")
        end
    end
    if mme.nModels > 1 && mme.M!=0
        if Pi != 0.0 && round(sum(values(Pi)),digits=2)!=1.0
          error("Summation of probabilities of Pi is not equal to one.")
        end
    end
end

function check_pedigree(mme,df,pedigree)
    if mme.ped == 0 && pedigree == false
        return
    elseif pedigree != false
        mme.ped = pedigree
    end
    if pedigree!=false
        pedID=map(string,collect(keys(pedigree.idMap)))
    else
        pedID=map(string,collect(keys(mme.ped.idMap)))
    end

    if mme.M!=0 && !issubset(mme.M.obsID,pedID)
        error("Not all genotyped individuals are found in pedigree!")
    end

    phenoID = map(string,df[1])
    if !issubset(phenoID,pedID)
        error("Not all phenotyped individuals are found in pedigree!")
    end
end

function check_outputID(mme)
    #Genotyped individuals are usaully not many, and are used in GWAS (complete
    #and incomplete), thus are used as default output_ID if not provided
    if mme.M == 0 && mme.pedTrmVec == 0
        mme.MCMCinfo.outputEBV = false #no EBV in non-genetic analysis
    end

    if mme.MCMCinfo.outputEBV == false
        mme.output_ID = 0
    elseif mme.output_ID == 0 && mme.M != 0
        mme.output_ID = mme.M.obsID
    elseif mme.output_ID == 0 && mme.M == 0 && mme.pedTrmVec != 0
        #output EBV for all individuals in the pedigree for PBLUP
        pedID=map(string,collect(keys(mme.ped.idMap)))
        mme.output_ID = pedID
    end
end

function check_phenotypes(mme,df)
    single_step_analysis = mme.MCMCinfo.single_step_analysis
    phenoID = map(string,df[1])   #same to df[:,1] in deprecated CSV
    if mme.M == 0 && mme.ped == 0 #non-genetic analysis
        return df
    end
    if single_step_analysis == false && mme.M != 0 #complete genomic data
        if !issubset(phenoID,mme.M.obsID)
            printstyled("Phenotyped individuals are not a subset of\n",
            "genotyped individuals (complete genomic data,non-single-step).\n",
            "Only use phenotype information for genotyped individuals.\n",bold=false,color=:red)
            index = [phenoID[i] in mme.M.obsID for i=1:length(phenoID)]
            df    = df[index,:]
            printstyled("The number of individuals with both genotypes and phenotypes used\n",
            "in the analysis is ",size(df,1),".\n",bold=false,color=:red)
        elseif mme.output_ID!=0 && !issubset(mme.output_ID,mme.M.obsID)
            printstyled("Testing individuals are not a subset of \n",
            "genotyped individuals (complete genomic data,non-single-step).\n",
            "Only output EBV for tesing individuals with genotypes.\n",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,mme.M.obsID)
        end
    else #incomplete genomic data , PBLUP
        pedID = map(string,collect(keys(mme.ped.idMap)))
        if !issubset(phenoID,pedID)
            printstyled("Phenotyped individuals are not a subset of\n",
            "individuals in pedigree (incomplete genomic data (single-step) or PBLUP).\n",
            "Only use phenotype information for individuals in the pedigree.\n",bold=false,color=:red)
            index = [phenoID[i] in pedID for i=1:length(phenoID)]
            df    = df[index,:]
            printstyled("The number of individuals with both phenotype and pedigree information\n",
            "used in the analysis is ",size(df,1),".\n",bold=false,color=:red)
        elseif mme.output_ID!=0 && !issubset(mme.output_ID,pedID)
            printstyled("Testing individuals are not a subset of \n",
            "individuals in pedigree (incomplete genomic data (single-step) or PBLUP).\n",
            "Only output EBV for tesing individuals in the pedigree.\n",bold=false,color=:red)
            mme.output_ID = intersect(mme.output_ID,pedID)
        end
    end
    return df
end

function init_mixed_model_equations(mme,df,sol)
    getMME(mme,df)
    #starting value for sol can be provided
    if sol == false #no starting values
        sol = zeros(size(mme.mmeLhs,1))
    else            #besure type is Float64
        nsol = mme.M!=0 ? size(mme.mmeLhs,1)+mme.M.nMarkers : size(mme.mmeLhs,1)
        if length(sol) != nsol || typeof(sol) <: AbstractString
            error("length or type of starting values is wrong.")
        end
        printstyled("Starting values are provided. The order of starting values for location parameters and\n",
        "marker effects should be the order of the Mixed Model Equation (This can be\n",
        "obtained by getNames(model)) and markers\n",bold=false,color=:red)
        sol = map(Float64,sol)
    end
    return sol,df
end
