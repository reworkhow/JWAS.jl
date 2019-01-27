"""
    add_genotypes(mme::MME,file,G;separator=' ',header=true,center=true,G_is_marker_variance=false,df=4.0)
* Get marker informtion from a genotype file.
* **G** defaults to the genetic variance with degree of freedom **df**=4.0.
* File format:

```
Animal,marker1,marker2,marker3,marker4,marker5
S1,1,0,1,1,1
D1,2,0,2,2,1
O1,1,2,0,1,0
O3,0,0,2,1,1
```
"""
function add_genotypes(mme::MME,file,G;separator=',',header=true,center=true,G_is_marker_variance=false,df=4)
    mme.M   = readgenotypes(file;separator=separator,header=header,center=center)
    if G_is_marker_variance == true
        mme.M.G = G
    else
        mme.M.genetic_variance = G
    end
    mme.df.marker = Float64(df)

    println(size(mme.M.genotypes,2), " markers on ",size(mme.M.genotypes,1)," individuals were added.")
end

"""
    same to add_genotypes
"""
function add_markers(mme::MME,file,G;separator=',',header=true,center=true,G_is_marker_variance=false,df=4)
    add_genotypes(mme,file,G;separator=separator,header=header,center=center,G_is_marker_variance=G_is_marker_variance,df=df)
end

#load genotypes from a text file (individual IDs in the 1st column)
function readgenotypes(file::AbstractString;separator=',',header=true,center=true)
    printstyled("The delimiter in ",split(file,['/','\\'])[end]," is \'",separator,"\'.\n",bold=false,color=:red)
    printstyled("The header (marker IDs) is ",(header ? "provided" : "not provided")," in ",split(file,['/','\\'])[end],".\n",bold=false,color=:red)

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
    genotypes = map(Float64,convert(Matrix,df[2:end]))
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
