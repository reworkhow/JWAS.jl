function make_yVecs(file::ASCIIString;header=false)
    df = readtable(file, eltypes=[UTF8String, Float64], separator = ' ',header=header)
    y  = df[:,2]
    return y
end

function make_yVecs(y::Array{Float64,1};header=false)
    y  = y
    return y
end