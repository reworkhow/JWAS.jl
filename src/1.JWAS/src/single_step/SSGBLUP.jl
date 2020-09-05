function get_Hi(pedigree,genotypes)
    #***********************************************************************
    #Hi = Ai + [0 0;0 inv(G)- inv(A_gg)]
    #***********************************************************************
    #reorder Ai to [non-genotyped ids; genotyped ids]
    num_pedn = PedModule.genoSet!(genotypes.obsID,pedigree)
    Ai       = PedModule.AInverse(pedigree)
    A        = inv(Array(Ai))
    A_gg     = A[(num_pedn+1):size(A,1),(num_pedn+1):size(A,1)]

    M,p,freq = genotypes.genotypes, genotypes.nMarkers, genotypes.alleleFreq #already centered
    M        = M ./ sqrt.(2*freq.*(1 .- freq))
    G        = (M*M'+ I*0.00001)/p

    IDs      = ped.IDs
    Hi       = Ai
    Hi[(num_pedn+1):size(Ai,1),(num_pedn+1):size(Ai,1)] += inv(G)-inv(A_gg)

    return Hi,IDs
end
