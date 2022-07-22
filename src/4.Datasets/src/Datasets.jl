module Datasets

using Printf,HTTP
function dataset(file_name::AbstractString;dataset_name::AbstractString="")
    basename = joinpath(dirname(@__FILE__), "..", "data", dataset_name)
    rdaname = joinpath(basename, string(file_name))
    if isfile(rdaname)
        return rdaname
    else
        error(@sprintf "Unable to locate file %s in %s\n" file_name dataset_name)
    end
end

function dataset_big(file_name::AbstractString;dataset_url::AbstractString="https://raw.githubusercontent.com/zhaotianjing/bio_protocol/main/data")
    rdaname = joinpath(dataset_url, string(file_name))
    http_obj = HTTP.get(rdaname)
    return http_obj
    end

export dataset


end # module
