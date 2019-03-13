#update
function getinfo(pedigree::PedModule.Pedigree)
    Ai   = PedModule.AInverse(pedigree)
    IDs = PedModule.getIDs(pedigree)
    inbreeding = PedModule.getInbreeding(pedigree)
    #outID=string.(1:12)
    #JWAS.mkmat_incidence_factor(outID,inID)
    #mat=Matrix(JWAS.mkmat_incidence_factor(outID,inID))
    #A=mat*A*mat'
    println("Get individual IDs, inverse of numerator relationship matrix, and inbreeding coefficients.")
    return IDs,Ai,inbreeding
end
