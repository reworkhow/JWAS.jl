module Datasets

using Printf
function dataset(file_name::AbstractString;dataset_name::AbstractString="")
    basename = joinpath(dirname(@__FILE__), "..", "data", dataset_name)
    rdaname = joinpath(basename, string(file_name))
    if isfile(rdaname)
        return rdaname
    else
        error(@sprintf "Unable to locate file %s in %s\n" file_name dataset_name)
    end
end

export dataset

end # module
