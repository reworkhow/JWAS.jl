"""
    add_genotypes(mme::MME,M::Union{AbstractString,Array{Float64,2},Array{Float32,2},Array{Any,2},DataFrames.DataFrame},G=false;
                  header=false,rowID=false,separator=',',
                  center=true,G_is_marker_variance=false,df=4)
* Get marker informtion from a genotype file or an nxp Matrix M of genotypes (Array or DataFrame),
  where n is the number of individuals and p is the number of markers. This file/matrix needs to be column-wise sorted by marker positions.
* **G** is the mean for the prior assigned for the genomic variance with degree of freedom **df**, defaulting to 4.0.
  If **G** is not provided, a value is calculated from responses (phenotypes)
* If a text file is provided, the file format should be:
  * ```
    Animal,marker1,marker2,marker3,marker4,marker5
    S1,1,0,1,1,1
    D1,2,0,2,2,1
    O1,1,2,0,1,0
    O3,0,0,2,1,1
    ```
* If an nxp Matrix of genotypes (Array or DataFrame) is provided, where n is the number of individuals and p is the number of markers,
  * This matrix needs to be column-wise sorted by marker positions.
  * rowID is a vector of individual IDs, e.g.,rowID=[\"a1\",\"b2\",\"c1\"]; if it is omitted, IDs will be set to 1:n
  * header is a header vector such as ["id"; "mrk1"; "mrk2";...;"mrkp"]. If omitted, marker names will be set to 1:p
"""
function add_genotypes(mme::MME,M::Union{AbstractString,Array{Float64,2},Array{Float32,2},Array{Any,2},DataFrames.DataFrame},G=false;
                       header=true,rowID=false,separator=',',
                       center=true,G_is_marker_variance=false,df=4)
    if G != false && size(G,1) != mme.nModels
        error("The covariance matrix is not a ",mme.nModels," by ",mme.nModels," matrix.")
    end
    genotypei   = get_genotypes(M,G,
                            header = header,rowID = rowID,separator = separator,
                            center = center,G_is_marker_variance = G_is_marker_variance,df = df)
    genotypei.ntraits = mme.nModels
    genotypei.trait_names = string.(mme.lhsVec)
    if mme.nModels != 1
      genotypei.df = genotypei.df + mme.nModels
    end
    genotypes = []
    push!(genotypes,genotypei)
    mme.M = genotypes
    if mme.nModels == 1 #?move to set_marker_hyperparameters_variances_and_pi?
        mme.df.marker = Float32(df)
    else
        νG0   = Float32(df) + mme.nModels
        mme.df.marker = νG0 #final df for inverse wishart
    end
end
################################################################################
#
# convert genotype file or matrix to struct Genotypes in MME
#
################################################################################
#1)load genotypes from a text file (1st column: individual IDs; 1st row: marker IDs (optional))
#2)load genotypes from Array or DataFrames (no individual IDs; no marker IDs (header))
"""
    get_genotypes(file::Union{AbstractString,Array{Float64,2},Array{Float32,2},Array{Any,2},DataFrames.DataFrame},G=false;
                  separator=',',header=true,rowID=false,
                  center=true,quality_control=false,
                  method = "RR-BLUP",Pi = 0.0,estimatePi = true,estimateScale=false,
                  G_is_marker_variance = false,df = 4.0)
* Get marker informtion from a genotype file/matrix. This file needs to be column-wise sorted by marker positions.
    * If a text file is provided, the file format should be:
      ```
      Animal,marker1,marker2,marker3,marker4,marker5
      S1,1,0,1,1,1
      D1,2,0,2,2,1
      O1,1,2,0,1,0
      O3,0,0,2,1,1
      ```
    * If an nxp Matrix of genotypes (Array or DataFrame) is provided, where n is the number of individuals and p is the number of markers,
        * This matrix needs to be column-wise sorted by marker positions.
        * rowID is a vector of individual IDs, e.g.,rowID=[\"a1\",\"b2\",\"c1\"]; if it is omitted, IDs will be set to 1:n
        * header is a header vector such as ["id"; "mrk1"; "mrk2";...;"mrkp"]. If omitted, marker names will be set to 1:p
* If `quality_control`=true, defaulting to `true`,
    * Missing genotypes should be denoted as `9`, and will be replaced by column means. Users can also impute missing genotypes before the analysis.
    * Minor allele frequency `MAF` threshold, defaulting to `0.01`, is uesd, and fixed loci are removed.
* **G** is the mean for the prior assigned for the genomic variance with degree of freedom **df**, defaulting to 4.0.
  If **G** is not provided, a value is calculated from responses (phenotypes).
* Available `methods` include "conventional (no markers)", "RR-BLUP", "BayesA", "BayesB", "BayesC", "Bayesian Lasso", and "GBLUP".
* In Bayesian variable selection methods, `Pi` for single-trait analyses is a number; `Pi` for multi-trait analyses is
  a dictionary such as `Pi=Dict([1.0; 1.0]=>0.7,[1.0; 0.0]=>0.1,[0.0; 1.0]=>0.1,[0.0; 0.0]=>0.1)`, defaulting to `all markers
  have effects (Pi = 0.0)` in single-trait analysis and `all markers have effects on all traits
  (Pi=Dict([1.0; 1.0]=>1.0,[0.0; 0.0]=>0.0))` in multi-trait analysis. `Pi` is estimated if `estimatePi` = true, , defaulting to `false`.
* Scale parameter for prior of marker effect variance is estimated if `estimateScale` = true, defaulting to `false`.

"""
function get_genotypes(file::Union{AbstractString,Array{Float64,2},Array{Float32,2},Array{Any,2},DataFrames.DataFrame},G=false;
                       method = "BayesC",Pi = 0.0,estimatePi = true, estimateVariance=true, estimateScale=false,
                       separator=',',header=true,rowID=false,
                       center=true,G_is_marker_variance = false,df = 4.0,
                       starting_value=false,
                       quality_control=true, MAF = 0.01, missing_value = 9.0)
    #Read the genotype file
    if typeof(file) <: AbstractString
        printstyled("The delimiterd in ",split(file,['/','\\'])[end]," is \'",separator,"\'. ",bold=false,color=:green)
        printstyled("The header (marker IDs) is ",(header ? "provided" : "not provided")," in ",split(file,['/','\\'])[end],".\n",bold=false,color=:green)
        #get marker IDs
        myfile = open(file)
        row1   = split(readline(myfile),[separator,'\n'],keepempty=false)
        if header==true
            markerID=string.(row1[2:end])  #skip header
        else
            markerID=string.(1:length(row1[2:end]))
        end
        #set type for each column
        ncol= length(row1)
        etv = Array{DataType}(undef,ncol)
        fill!(etv,Float32)
        etv[1]=String
        close(myfile)
        #read a large genotype file
        data      = CSV.read(file,DataFrame,types=etv,delim = separator,header=false,skipto=(header==true ? 2 : 1))
        obsID     = data[!,1]
        genotypes = Matrix(data[!,2:end])
    elseif typeof(file) == DataFrames.DataFrame #Datafarme
        println("The first column in the dataframe should be individual IDs.")
        println("The data type of markers should be Number.")
        if header == true
            markerID = names(file)[2:end]
        else
            markerID = string.(1:(size(file,2)-1))
        end
        obsID     = map(string,file[!,1])
        genotypes = map(Float32,Matrix(file[!,2:end]))
    elseif typeof(file) <: Union{Array{Float64,2},Array{Float32,2},Array{Any,2}} #Array
        if length(header) != (size(file,2)+1)
            header = ["id"; string.(1:size(file,2))]
            printstyled("The marker IDs are set to 1,2,...,#markers\n",bold=true)
        end
        if length(rowID) != size(file,1)
            rowID = string.(1:size(file,1))
            printstyled("The individual IDs is set to 1,2,...,#observations\n",bold=true)
        end
        markerID  = string.(header[2:end])
        obsID     = map(string,rowID)
        genotypes = map(Float32,convert(Matrix,file))
    else
        error("The data type is not supported.")
    end

    isGRM = false
    if method == "GBLUP" #fix and check some parameters in GBLUP
        if issymmetric(genotypes) #a kernel / relationship matrix is provided
            center = false
            quality_control = false
            isGRM = true
            println("A genomic relationship matrix is provided (instead of a genotype covariate matrix).")
        end
        if G_is_marker_variance == true
            error("Genetic variance is required.")
        end
        if starting_value != false
            error("starting values are not supported for GBLUP now.")
        end
    end

    #preliminary summary of genotype
    nObs,nMarkers = size(genotypes)       #number of individuals and molecular markers

    #Naive Quality Control 1, replace missing values with column means
    if quality_control == true
        for genoi in eachcol(genotypes)
            missing_obs        = findall(x->x==float(missing_value),genoi)
            nonmissing_obs     = deleteat!(collect(1:nObs),missing_obs)
            genoi[missing_obs] .= mean(genoi[nonmissing_obs])
            if findfirst(x->(x>2.0||x<0.0),genoi) != nothing #issue71, genotype score
                @warn "genotype scores out of the range 0 to 2 are found."
            end
        end
        printstyled("Missing values ($missing_value) are replaced by column means.\n",bold=true)
    end

    markerMeans   = center==true ? center!(genotypes) : mean(genotypes,dims=1) #centering genotypes or not
    p             = markerMeans/2.0       #allele frequency

    #Naive Quality Control 2, minor allel frequency & fixed loci
    if quality_control == true
         select1   = MAF .< vec(p) .< 1-MAF
         select2   = vec(var(genotypes,dims=1)) .!= 0
         select    = select1 .& select2
         genotypes = genotypes[:,select]
         p         = p[:,select]
         markerID  = markerID[select]
         printstyled("$(sum(1 .- select)) loci which are fixed or have minor allele frequency < $MAF are removed.\n",bold=true)
    end
    nObs,nMarkers = size(genotypes)       #number of individuals and molecular markers
    sum2pq        = (2*p*(1 .- p)')[1,1]  #∑2pq

    #a kernel / relationship matrix or genotypes are provided in GBLUP
    if method == "GBLUP"
        if !isGRM #calculate the relationship matrix from the genotype covariate matrix
            genotypes  = genotypes ./ sqrt.(2*p.*(1 .- p))
            genotypes  = (genotypes*genotypes'+ I*0.00001)/nMarkers
            println("A genomic relationship matrix is computed from genotypes.")
            isGRM  = true
        end
        add_small_value_count = 0
        while isposdef(genotypes) == false
            println("The relationship matrix is not positive definite. A very small number is added to the diagonal.")
            genotypes = genotypes + I*0.00001
            add_small_value_count += 1
            if add_small_value_count > 10
                error("Please provide a positive-definite relationship matrix.")
            end
        end
    end

    genotypes     = Genotypes(obsID,markerID,nObs,nMarkers,p,sum2pq,center,genotypes,isGRM)
    if G_is_marker_variance == true
        genotypes.G = G
    else
        genotypes.genetic_variance = G
    end
    genotypes.method     = method
    genotypes.estimatePi = estimatePi
    genotypes.π          = Pi
    genotypes.df         = df #It will be modified base on number of traits in build_model()
    genotypes.estimateScale    = estimateScale
    genotypes.estimateVariance = estimateVariance

    writedlm("IDs_for_individuals_with_genotypes.txt",genotypes.obsID)
    println("Genotype informatin:")
    println("#markers: ",(isGRM ? 0 : size(genotypes.genotypes,2)),"; #individuals: ",size(genotypes.genotypes,1))

    #starting values for marker effects
    if starting_value != false
        printstyled("Starting values are provided for marker efffects. The order of starting values should be\n",
        "markers for each traits (all markers for trait 1 then all markers for trait 2...)\n",bold=false,color=:green)
        nsol = (method != "GBLUP" ? nMarkers : nObs)
        if length(starting_value)%nsol != 0 #work for both single and multiple traits
            error("length of starting values is wrong.")
        end
        genotypes.α = starting_value
    end
    return genotypes
end
