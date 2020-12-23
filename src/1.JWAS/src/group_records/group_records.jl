# ID,group,y1,y2,y3,x1,x2,x3,dam
# a1,group1,-0.06,3.58,-1.18,0.9,2,m,a2
# a3,group2,-2.07,3.19,1,0.7,2,f,a3
# a4,group1,-2.63,6.97,-0.83,0.6,1,m,a2
# a5,group1,2.31,5.3,-1.52,0.4,2,m,a2
# a6,group2,1.5,3,2,0.5,1,f,a3
# a7,group1,2,4.5,0.8,0.1,0.1,m,a3
# a8,group2,-1.4,1,1,0.3,1,f,a6

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
