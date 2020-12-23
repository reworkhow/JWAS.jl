# ID,group_records,y1
# a1,group1,-0.06
# a3,group2,-2.07
# a4,group1,-2.63
# a5,group1,2.31
# a6,group2,1.5
# a7,group1,2.0
# a8,group2,-1.4

using StatsBase,Statistics,Printf,Random,DelimitedFiles,InteractiveUtils,DataFrames,CSV,SparseArrays,LinearAlgebra
phenotypes = CSV.read("/Users/lijinghui/Desktop/Genetics/Group_GWAS/function_test/phenotypes.csv",
                delim = ',',header=true,missingstrings=["NA"])

#the column name for group records hs to be "group_records"
function get_T(data)
    uID = data[:,"group_records"] # group info
    yID = unique(uID) #unique groups

    #get an incidence matrix Z to reorder uID to yID by yID = Z*uID
    T = mkmat_incidence_factor(yID,uID)
    return T
end
