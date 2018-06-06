################################################################################
#User-interface to output MCMC Samples for specific variables
################################################################################
"""
    outputEBV(model,IDs::Array{AbstractString,1})

Output estimated breeding values and prediction error variances for IDs.
"""
function outputEBV(model,IDs::Array{AbstractString,1})
    println("Estimated breeding values and prediction error variances will be included in the output.")
end
