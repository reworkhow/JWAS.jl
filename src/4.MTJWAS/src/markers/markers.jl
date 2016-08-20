include("tools.jl")
include("readgenotypes.jl")
include("BayesianAlphabet/BayesC0.jl")
include("BayesianAlphabet/BayesC.jl")
include("BayesianAlphabet/BayesB.jl")
include("BayesianAlphabet/MTBayesC0.jl")
include("BayesianAlphabet/MTBayesC.jl")
include("BayesianAlphabet/MTBayesCC.jl")
include("BayesianAlphabet/MTBayesB.jl")

function addMarkers(mme::MME,file,G;
                    separator=' ',header=true,
                    G_is_marker_variance=false)
    mme.M   = readgenotypes(file;separator=separator,header=header,center=true)
    mme.M.G = G
    mme.M.G_is_marker_variance = G_is_marker_variance
end
