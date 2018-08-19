function readgenotypes(file::AbstractString;separator=' ',header=false,center=true)
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


    df            = readdlm(file,separator,header=header)

    if header == true
        df=df[1]
    end
    obsID         = map(String,df[:,1])
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

#M is of type: Array{Float64,2}
function readgenotypes(M::Union{Array{Float64,2},DataFrames.DataFrame};header=false,separator=' ',center=true)
    header        = false
    obsID         = map(String,M[:,1])
    markerID      = ["NA"]
    genotypes     = map(Float64,M[:,2:end])
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
