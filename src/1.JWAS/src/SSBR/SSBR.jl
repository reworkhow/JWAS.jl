function SSBRrun(mme,df)
    obsID    = map(string,df[!,1])          #phenotyped ID
    geno     = mme.M                        #input genotyps
    ped      = mme.ped                      #pedigree
    Ai_nn,Ai_ng = calc_Ai(ped,geno,mme)     #get A inverse
    impute_genotypes(geno,ped,mme,Ai_nn,Ai_ng) #impute genotypes for non-genotyped inds

    #add model terms for SSBR
    add_term(mme,"ϵ") #impuatation residual
    add_term(mme,"J") #centering parameter
    #add data for ϵ and J (add columns in input phenotypic data)
    isnongeno = [ID in mme.ped.setNG for ID in obsID] #true/false
    data_ϵ    = deepcopy(map(string,df[!,1]))
    data_ϵ[.!isnongeno].="0"
    df[!,Symbol("ϵ")]=data_ϵ

    df[!,Symbol("J")],mme.output_X["J"]=make_JVecs(mme,df,Ai_nn,Ai_ng)
    set_covariate(mme,"J")
    set_random(mme,"ϵ",mme.M.genetic_variance,Vinv=Ai_nn,names=ped.IDs[1:size(Ai_nn,1)])
    if mme.M.genetic_variance == false
        error("Please input the genetic variance using add_genotypes()")
    end

    outputMCMCsamples(mme,"J")
    outputMCMCsamples(mme,"ϵ")
end

############################################################################
# Phenotypes (any order )
############################################################################

############################################################################
# Pedigree
############################################################################
function calc_Ai(ped,geno,mme)
    #***********************************************************************
    #***reorder in A (ped) as genotyped id then non-genotyped id
    #***********************************************************************
    num_pedn = PedModule.genoSet!(geno.obsID,ped)
    #calculate Ai, Ai_nn, Ai_ng
    mme.Ai   = PedModule.AInverse(ped)
    Ai_nn    = mme.Ai[1:num_pedn,1:num_pedn]
    Ai_ng    = mme.Ai[1:num_pedn,(num_pedn+1):size(mme.Ai,1)]
    return Ai_nn,Ai_ng
end
############################################################################
# Genotypes
############################################################################
function impute_genotypes(geno,ped,mme,Ai_nn,Ai_ng)
    num_pedn = size(Ai_nn,1)
    #reorder genotypes to get Mg with the same order as Ai_gg
    Z  = mkmat_incidence_factor(ped.IDs[(num_pedn+1):end],geno.obsID)
    Mg = Z*geno.genotypes
    #impute genotypes for non-genotyped inds
    Mfull = [Ai_nn\(-Ai_ng*Mg);Mg]; #May use big memory, solve it using chunks of Mg

    mme.M.genotypes = Mfull
    mme.M.obsID     = ped.IDs
    mme.M.nObs      = length(mme.M.obsID)
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
        Zo  = mkmat_incidence_factor(mme.output_ID,mme.M.obsID)
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
