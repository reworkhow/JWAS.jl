#load genotypes from a text file (individual IDs in the 1st column)
function readgenotypes(file::AbstractString;separator=',',header=true,center=true)
    printstyled("The delimiter in ",file," is ",bold=false)
    printstyled("\'",separator,"\'.\n",bold=false,color=:red)


    myfile = open(file)
    #get number of columns
    row1   = split(readline(myfile),[separator,'\n'],keepempty=false)

    if header==true
      markerID=row1[2:end]  #skip header
    else
      markerID= ["NA"]
    end
    #set types for each column and get number of markers
    ncol= length(row1)
    etv = Array{DataType}(undef,ncol)
    fill!(etv,Float64)
    etv[1]=String
    close(myfile)
    #
    # #read genotypes
    #df = CSV.read(file, types=etv, delim = separator, header=header)
    df = readtable(file, eltypes=etv, separator = separator, header=header)
    obsID     = map(String,df[1]) #convert from Array{Union{String, Missings.Missing},1} to String #redundant actually
    genotypes = map(Float64,convert(Array,df[2:end]))
    nObs,nMarkers = size(genotypes)

    ##readdlm
    #df            = readdlm(file,separator,Any,header=header)
    #if header == true
    #    df=df[1]
    #end
    #obsID         = map(string,df[:,1]) #map(String,df[:,1]) not work, if df[:,1] are integers
    #genotypes     = map(Float64,df[:,2:end])
    #nObs,nMarkers = size(genotypes)

    if center==true
        markerMeans = center!(genotypes) #centering
    else
        markerMeans = center!(copy(genotypes))  #get marker means
    end
    p             = markerMeans/2.0
    sum2pq        = (2*p*(1 .- p)')[1,1] #(1x1 matrix)[1,1]

    return Genotypes(obsID,markerID,nObs,nMarkers,p,sum2pq,center,genotypes)
end

#load genotypes from Array or DataFrames (individual IDs in the 1st column,
function readgenotypes(df::Union{Array{Float64,2},DataFrames.DataFrame};header=false,separator=false,center=true)
    if header==true
        error("The header (marker IDs) must be false or provided as an array when the input is not a file.")
    end

    if length(header)==size(df,2)
        markerID = header
    elseif header==false
        markerID = ["NA"]
    else
        error("The length of the header (marker IDs) must be equal to the number of markers (columns).")
    end

    if typeof(df) == DataFrames.DataFrame
        obsID     = map(string,df[1])
        genotypes = map(Float64,convert(Array,df[2:end]))
    else
        obsID         = map(string,view(df,:,1))  #map(String,df[:,1]) not work, if df[:,1] are integers
        genotypes     = map(Float64,convert(Array,view(df,:,2:size(df,2))))
    end
    nObs,nMarkers = size(genotypes)

    if center==true
        markerMeans = center!(genotypes) #centering
    else
        markerMeans = center!(copy(genotypes))  #get marker means
    end
    p             = markerMeans/2.0
    sum2pq       = (2*p*(1 .- p)')[1,1]

    return Genotypes(obsID,markerID,nObs,nMarkers,p,sum2pq,center,genotypes)
end
