include("tools4genotypes.jl")
include("readgenotypes.jl")
include("BayesianAlphabet/BayesC0.jl")
include("BayesianAlphabet/BayesC.jl")
include("BayesianAlphabet/BayesB.jl")
include("BayesianAlphabet/MTBayesC.jl")
include("BayesianAlphabet/MTBayesCC.jl")
include("BayesianAlphabet/MTBayesB.jl")
include("BayesianAlphabet/MTBayesC0L.jl")


"""
    add_genotypes(mme::MME,file,G;separator=' ',header=true,center=true,G_is_marker_variance=false,df=4.0)
* Get marker informtion from a genotype file (same order as the phenotype file).
* **G** defaults to the genetic variance with degree of freedom **df**=4.0.
* File format:

```
Animal,marker1,marker2,marker3,marker4,marker5
S1,1,0,1,1,1
D1,2,0,2,2,1
O1,1,2,0,1,0
O3,0,0,2,1,1
```
"""
function add_genotypes(mme::MME,file,G;separator=' ',header=true,center=true,rowID=true,G_is_marker_variance=false,df=4)
    mme.M   = readgenotypes(file;separator=separator,header=header,rowID=rowID,center=center)
    mme.M.G = G
    mme.M.G_is_marker_variance = G_is_marker_variance
    mme.df.marker = Float64(df)

    println(size(mme.M.genotypes,2), " markers on ",size(mme.M.genotypes,1)," individuals were added.")
end

"""
    same to add_genotypes
"""
function add_markers(mme::MME,file,G;separator=' ',header=true,center=true,rowID=true,G_is_marker_variance=false,df=4)
    add_genotypes(mme,file,G;separator=separator,header=header,center=center,G_is_marker_variance=G_is_marker_variance,df=df)
end
