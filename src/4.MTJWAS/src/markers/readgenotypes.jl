function readgenotypes(file::AbstractString;separator=' ',header=false,center=true)
    println("The delimiters in file $file is ",separator,"  .")

    myfile = open(file)
    #get number of columns
    row1   = split(readline(myfile),[separator,'\n'],keep=false)

    if header==true
      markerID=row1[2:end]  #skip header
    else
      markerID= ["NA"]
    end
    #set types for each column and get number of markers
    ncol= length(row1)
    etv = Array(DataType,ncol)
    fill!(etv,Float64)
    etv[1]=UTF8String
    close(myfile)

    #read genotypes
    df = readtable(file, eltypes=etv, separator = separator, header=header)
    #df = readtable(file, separator = separator, header=header)

    obsID=convert(Array,df[:,1])
    genotypes = map(Float64,convert(Array,df[:,2:end]))
    nObs,nMarkers = size(genotypes)

    if center==true
        markerMeans = center!(genotypes) #centering
    else
        markerMeans = center!(copy(genotypes))  #get marker means
    end
    p             = markerMeans/2.0
    sum2pq        = (2*p*(1-p)')[1,1]

    return Genotypes(obsID,markerID,nObs,nMarkers,p,sum2pq,center,genotypes)
end

#M is of type: Array{Float64,2} or DataFrames
function readgenotypes(M::DataFrames.DataFrame;header=false,separator=' ',center=true)#improve later
    genotypes  = map(Float64,convert(Array,M))
    obsID     = ["NA"]
    markerID  = ["NA"]
    nObs,nMarkers = size(genotypes)

    if center==true
        markerMeans = center!(genotypes) #centering
    else
        markerMeans = center!(copy(genotypes))  #get marker means
    end
    p             = markerMeans/2.0
    sum2pq       = (2*p*(1-p)')[1,1]

    return Genotypes(obsID,markerID,nObs,nMarkers,p,sum2pq,center,genotypes)
end

function readgenotypes(M::Array{Float64,2};header=false,separator=' ',center=true)
    genotypes = map(Float64,M)
    obsID     = ["NA"]
    markerID  = ["NA"]
    nObs,nMarkers = size(genotypes)

    if center==true
        markerMeans = center!(genotypes) #centering
    else
        markerMeans = center!(copy(genotypes))  #get marker means
    end
    p             = markerMeans/2.0
    sum2pq       = (2*p*(1-p)')[1,1]

    return Genotypes(obsID,markerID,nObs,nMarkers,p,sum2pq,center,genotypes)
end

export readgenotypes
