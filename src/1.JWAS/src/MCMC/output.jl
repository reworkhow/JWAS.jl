################################################################################
#User-interface to output EBV
################################################################################
"""
    outputEBV(model,IDs::Array{String,1})

Output estimated breeding values and prediction error variances for IDs.
"""
function outputEBV(model,IDs::Array{String,1})
    println("Estimated breeding values and prediction error variances will be included in the output.")
    model.output_ID=IDs
end

function reformat(res::Array)
    out_names=[strip(i) for i in split(res[1,1],':',keep=false)]
    for rowi in 2:size(res,1)
        out_names=[out_names [strip(i) for i in split(res[rowi,1],':',keep=false)]]
    end
    out_names=permutedims(out_names,[2,1])
    out_values=map(Float64,res[:,2])
    out=[out_names out_values]
    out = DataFrame(out, [:Trait, :Effect, :Level, :Estimate],)
    return out
end
