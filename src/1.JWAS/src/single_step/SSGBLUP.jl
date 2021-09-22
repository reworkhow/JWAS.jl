function get_Hi(pedigree,genotypes)
    #***********************************************************************
    #Hi = Ai + [0 0;0 inv(G)- inv(A_gg)]
    #***********************************************************************
    #reorder Ai to [non-genotyped ids; genotyped ids]
    num_pedn = PedModule.genoSet!(genotypes.obsID,pedigree)
    Ai       = PedModule.AInverse(pedigree)
    A        = inv(Array(Ai))
    A_gg     = A[(num_pedn+1):size(A,1),(num_pedn+1):size(A,1)]


    if genotypes.method == "GBLUP"
        G = genotypes.genotypes
    else
        M,p,freq = genotypes.genotypes, genotypes.nMarkers, genotypes.alleleFreq #already centered
        M        = M ./ sqrt.(2*freq.*(1 .- freq))
        G        = (M*M'+ I*0.00001)/p
        add_small_value_count = 0
        while isposdef(G) == false
            println("The relationship matrix is not positive definite. A very small number is added to the diagonal.")
            G += I*0.00001
            add_small_value_count += 1
            if add_small_value_count > 10
                error("Please provide a positive-definite realtionship matrix.")
            end
        end
    end

    IDs      = pedigree.IDs
    Hi       = Ai
    Hi[(num_pedn+1):size(Ai,1),(num_pedn+1):size(Ai,1)] += inv(G)-inv(A_gg)

    return Hi,IDs
end
