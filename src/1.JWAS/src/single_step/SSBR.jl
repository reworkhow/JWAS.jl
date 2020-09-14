function SSBRrun(mme,df,big_memory=false)
    obsID    = mme.obsID                    #phenotyped ID
    geno     = mme.M[1]                     #input genotyps
    ped      = mme.ped                      #pedigree
    println("calculating A inverse")
    flush(stdout)
    @time Ai_nn,Ai_ng = calc_Ai(ped,geno,mme)     #get A inverse
    println("imputing missing genotypes")
    flush(stdout)
    @time impute_genotypes(geno,ped,mme,Ai_nn,Ai_ng,big_memory) #impute genotypes for non-genotyped inds
    println("completed imputing genotypes")
    #add model terms for SSBR
    add_term(mme,"ϵ") #impuatation residual
    if mme.MCMCinfo.fitting_J_vector == true
        add_term(mme,"J") #centering parameter
    end
    #add data for ϵ and J (add columns in input phenotypic data)
    isnongeno = [ID in mme.ped.setNG for ID in obsID] #true/false
    data_ϵ    = deepcopy(obsID)
    data_ϵ[.!isnongeno].="missing"
    df[!,Symbol("ϵ")]=data_ϵ

    if mme.MCMCinfo.fitting_J_vector == true
        df[!,Symbol("J")],mme.output_X["J"]=make_JVecs(mme,df,Ai_nn,Ai_ng)
        set_covariate(mme,"J")
    end
    set_random(mme,"ϵ",geno.genetic_variance,Vinv=Ai_nn,names=ped.IDs[1:size(Ai_nn,1)])
    if geno.genetic_variance == false
        error("Please input the genetic variance using add_genotypes()")
    end
    if mme.MCMCinfo.fitting_J_vector == true
        outputMCMCsamples(mme,"J")
    end
    outputMCMCsamples(mme,"ϵ")
end

############################################################################
# Phenotypes (any order )
############################################################################

############################################################################
# Pedigree
############################################################################
function calc_Ai(ped,geno,mme)
    #reorder Ai to [non-genotyped ids; genotyped ids]
    num_pedn = PedModule.genoSet!(geno.obsID,ped)
    #calculate Ai, Ai_nn, Ai_ng
    Ai       = PedModule.AInverse(ped)
    Ai_nn    = Ai[1:num_pedn,1:num_pedn]
    Ai_ng    = Ai[1:num_pedn,(num_pedn+1):size(Ai,1)]
    return Ai_nn,Ai_ng
end
############################################################################
# Genotypes
############################################################################
function impute_genotypes(geno,ped,mme,Ai_nn,Ai_ng,big_memory=false)
    num_pedn = size(Ai_nn,1)
    #reorder genotypes to get Mg with the same order as Ai_gg
    Z  = mkmat_incidence_factor(ped.IDs[(num_pedn+1):end],geno.obsID)
    Mg = Z*geno.genotypes
    #**********************************************************
    #impute genotypes for non-genotyped inds
    #**********************************************************
    #two approached are available with similar speed and memory
    #appraoch 1
    #Ai_nn*Mn = -Ai_ng*Mg => [Ai_nn\(-Ai_ng*Mg)
    #                         Mg];
    #
    #appraoch 2
    #[Ai_nn 0 [Mn    [-Ai_ng
    #    0 I]* Mg] =       I]*Mg
    ngeno   = size(Ai_ng,2)
    block12 = spzeros(num_pedn, ngeno)
    block21 = spzeros(ngeno,num_pedn)
    block22 = sparse(I, ngeno, ngeno)
    lhs     = [Ai_nn block12
               block21 block22];
    rhs     = [-Ai_ng
               block22];
    #genotypes are imputed for all individuals in pedigree and aligned to
    #phenotyped individuals by genotype chunks to avoid memory issues.
    Z              = mkmat_incidence_factor(mme.obsID,ped.IDs)  #aligned to phenotypes
    Zo             = mkmat_incidence_factor(mme.output_ID,ped.IDs) #aligned to outputIDs
    if big_memory == false
        nmarkers       = size(Mg,2)
        Mpheno         = zeros(length(mme.obsID),nmarkers)
        Mout           = zeros(length(mme.output_ID),nmarkers)
        markerperchunk = 1000
        if nmarkers%markerperchunk == 0
            nchunks = div(nmarkers,markerperchunk)
        else
            nchunks = div(nmarkers,markerperchunk)+1
        end
        @showprogress "imputing genotypes for non-genotyped individuals" for i=1:nchunks
            if i != nchunks
                myrange = (1+(i-1)*markerperchunk):(i*markerperchunk)
            else
                myrange = (1+(i-1)*markerperchunk):nmarkers
            end
            Mped_chunk        = lhs\(rhs*Mg[:,myrange])
            Mpheno[:,myrange] = Z*Mped_chunk
            Mout[:,myrange]   = Zo*Mped_chunk
        end
    else
         Mpheno = Z*(lhs\(rhs*Mg))
         Mout   = Zo*(lhs\(rhs*Mg))
    end

    geno.output_genotypes = Mout
    geno.genotypes        = Mpheno
    geno.obsID            = mme.obsID
    geno.nObs             = length(mme.obsID)
    GC.gc()
end
############################################################################
# Fixed effects (J)
############################################################################
function make_JVecs(mme,df,Ai_nn,Ai_ng)
    mme.obsID =  map(string,df[!,1])
    Jg = -ones(size(Ai_ng,2))
    Jn = Ai_nn\(-Ai_ng*Jg)
    J  = [Jn;
          Jg]
    Z  = mkmat_incidence_factor(mme.obsID,mme.ped.IDs)
    if mme.output_ID != 0
        Zo  = mkmat_incidence_factor(mme.output_ID,mme.ped.IDs)
        ZoJ = Zo*J
    else
        ZoJ=0
    end
    return Z*J,ZoJ
end

"""
add to model an extra term: imputation_residual
"""
function add_term(mme,term::AbstractString)
    for m in 1:mme.nModels
        modelterm = ModelTerm(term,m,string(mme.lhsVec[m]))
        push!(mme.modelTerms,modelterm)
        mme.modelTermDict[modelterm.trmStr]=modelterm
    end
end
