#http://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=singlestepblupf90.pdf
function get_Hi(pedigree,genotypes;weight_for_G::Float64 = 1.0)
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
    end
    if 0.0 <= weight_for_G <= 1.0
        weight_for_G = round(weight_for_G, digits=3)
        weight_for_A = round(1-weight_for_G, digits=3)
        G = weight_for_G*G + weight_for_A*A_gg
        println("The genomic relationship matrix is calculated as
        G = $weight_for_G G + $weight_for_A A.")
    end
    if isposdef(G) == false
        error("Please provide a positive-definite relationship matrix.")
    end

    IDs      = pedigree.IDs
    Hi       = Ai
    Hi[(num_pedn+1):size(Ai,1),(num_pedn+1):size(Ai,1)] += inv(G)-inv(A_gg)

    return Hi,IDs
end
