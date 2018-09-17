function readgenotypes(file::AbstractString;separator=' ',header=false,rowID=true,center=true)
    #println("The delimiters in file $file is ",separator,"  .")

    myfile = open(file)
    #get number of columns
    row1   = split(readline(myfile),[separator,'\n'],keepempty=false)

    if header==true
      markerID=row1[2:end]  #skip header
    else
      markerID= ["NA"]
    end
    #set types for each column and get number of markers
    # ncol= length(row1)
    # etv = Array{DataType}(ncol)
    # fill!(etv,Float64)
    # etv[1]=String
    # close(myfile)
    #
    # #read genotypes
    # df = CSV.read(file, types=etv, delim = separator, header=header)
    # obsID     = map(String,df[:,1]) #convert from Array{Union{String, Missings.Missing},1} to String #redundant actually
    # genotypes = map(Float64,convert(Array,df[:,2:end]))
    # nObs,nMarkers = size(genotypes)


    df            = readdlm(file,separator,Any,header=header)

    if header == true
        df=df[1]
    end
    obsID         = map(string,df[:,1]) #map(String,df[:,1]) not work, if df[:,1] are integers
    genotypes     = map(Float64,df[:,2:end])
    nObs,nMarkers = size(genotypes)



    if center==true
        markerMeans = center!(genotypes) #centering
    else
        markerMeans = center!(copy(genotypes))  #get marker means
    end
    p             = markerMeans/2.0
    sum2pq        = (2*p*(1 .- p)')[1,1]

    return Genotypes(obsID,markerID,nObs,nMarkers,p,sum2pq,center,genotypes)
end

#M is of type: Array{Float64,2} or DataFrames ( genotype covariates only), marker ID and row ID required separately.
function readgenotypes(M::Union{Array{Float64,2},DataFrames.DataFrame};header=true,rowID=true,separator=' ',center=true)
    if rowID == true
        @error "Row IDs must be provided as an array when the input is not a file."
    end
    if header==true
        @error "Header (marker IDs) must be false or provided as an array when the input is not a file."
    end

    if length(rowID)==size(M,1)
        obsID      = rowID
    else
        @error "The length of row IDs must be equal to the number of individuals in the genotype covariate matrix."
    end
    if length(header)==size(M,2)
        markerID = header
    elseif header==false
        markerID = ["NA"]
    else
        @error "The length of header (marker IDs) must be equal to the number of markers in the genotype covariate matrix."
    end

    genotypes  = map(Float64,convert(Array,M))


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
