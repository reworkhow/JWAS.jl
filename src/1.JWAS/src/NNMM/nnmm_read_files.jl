# function read_genotypes(file::Union{AbstractString,Array{Float64,2},Array{Float32,2},Array{Int64,2}, Array{Int32,2}, Array{Any,2}, DataFrames.DataFrame}, G = false;
#                     #    ## method:
#                     #    method = "BayesC", Pi = 0.0, estimatePi = true, 
#                     #    ## variance:
#                     #    G_is_marker_variance = false, df = 4.0,
#                     #    estimate_variance = true, estimate_scale = false,
#                     #    constraint = false, #for multi-trait only, constraint=true means no genetic covariance among traits
#                        ## format:
#                        separator = ',', header = true, double_precision = false,
#                        ## quality control:
#                        quality_control = true, MAF = 0.01, missing_value = 9.0,
#                        ## others:
#                        center = true, starting_value = false)
#     #Define data precision
#     data_type = double_precision ? Float64 : Float32 #works for both genotypes and GRM
#     #Read the genotype file
#     if typeof(file) <: AbstractString  #string (path to a file)
#         #print info about the file format
#         printstyled("The delimiter in ", split(file, ['/', '\\'])[end], " is \'", separator, "\'. ", bold=false, color=:green)
#         printstyled("The header (marker IDs) is ", (header ? "provided" : "not provided"), " in ", split(file, ['/', '\\'])[end], ".\n", bold=false, color=:green)
#         #get marker IDs
#         myfile = open(file)
#         row1   = split(readline(myfile), [separator, '\n'], keepempty=false) #read the first row
#         if header == true
#             markerID = string.(row1[2:end])           #extracts marker IDs from row1
#         else
#             markerID = string.(1:length(row1[2:end])) #generate default marker IDs if no header
#         end
#         #set type for each column
#         ncol = length(row1)
#         etv  = Array{DataType}(undef, ncol)
#         fill!(etv, data_type)
#         etv[1] = String #individual ID
#         close(myfile)
#         #read a large genotype file
#         data      = CSV.read(file,DataFrame,types=etv,delim = separator,header=false,skipto=(header==true ? 2 : 1))
#         obsID     = map(string,data[!,1])
#         genotypes = Matrix(data[!,2:end]) #memory usage is doubled here
#         #clean memory
#         data = nothing
#         GC.gc()
#         #Note: why we choose to use `genotypes = Matrix(data[!,2:end])`
#          #- cannot use `genotypes=@view data[!,2:end]` because the view will be copied as a new DataFrame in QC: `genotypes = genotypes[:,select]`
#          #- cannot use `select!(genotypes, Not(1))` to simply delete the 1st column of obsID from the DataFrame, because genotypes must be a Matrix used for matrix multiplication (e.g. EBV += Mi.output_genotypes*Mi.α[traiti]).
#     elseif typeof(file) == DataFrames.DataFrame #Datafarme
#         println("The first column in the dataframe should be individual IDs.")
#         println("The remaining columns are markers with the data type Number.")
#         if header == true
#             markerID = names(file)[2:end]
#         else
#             markerID = string.(1:(size(file,2)-1))
#             printstyled("The marker IDs are set to 1,2,...,#markers\n",bold=true)
#         end
#         obsID     = map(string,file[!,1])
#         genotypes = map(data_type, Matrix(file[!,2:end]))
#     elseif typeof(file) <: Union{Array{Float64,2}, Array{Float32,2}, Array{Int64,2}, Array{Int32,2}, Array{Any,2}} #Array (Matrix)
#         println("The input data is a genotype matrix, without individual IDs and marker IDs.")
#         markerID  = string.(1:size(file,2))
#         printstyled("The marker IDs are set to 1,2,...,#markers\n",bold=true)
#         obsID     = map(string,string.(1:size(file,1)))
#         printstyled("The individual IDs is set to 1,2,...,#observations\n",bold=true)
#         genotypes = map(data_type, convert(Matrix,file))
#     else
#         error("The data type is not supported.")
#     end

#     isGRM = false

#     # if method == "GBLUP" #fix and check some parameters in GBLUP
#     #     if issymmetric(genotypes) #a kernel / relationship matrix is provided
#     #         center = false
#     #         quality_control = false
#     #         isGRM = true
#     #         println("A genomic relationship matrix is provided (instead of a genotype covariate matrix).")
#     #     end
#     #     if G_is_marker_variance == true
#     #         error("Genetic variance is required.")
#     #     end
#     #     if starting_value != false
#     #         error("starting values are not supported for GBLUP now.")
#     #     end
#     # end

#     #preliminary summary of genotype
#     nObs,nMarkers = size(genotypes)       #number of individuals and molecular markers

#     #Naive Quality Control 1, replace missing values with column means
#     if quality_control == true
#         for genoi in eachcol(genotypes)
#             missing_obs        = findall(x->x==float(missing_value),genoi)
#             nonmissing_obs     = deleteat!(collect(1:nObs),missing_obs)
#             genoi[missing_obs] .= mean(genoi[nonmissing_obs])
#             if findfirst(x->(x>2.0||x<0.0),genoi) != nothing #issue71, genotype score
#                 @warn "genotype scores out of the range 0 to 2 are found."
#             end
#         end
#         printstyled("Missing values ($missing_value) are replaced by column means.\n",bold=true)
#     end

#     markerMeans   = center==true ? center!(genotypes) : mean(genotypes,dims=1) #centering genotypes or not
#     p             = markerMeans/data_type(2.0)       #allele frequency

#     #Naive Quality Control 2, minor allel frequency & fixed loci
#     if quality_control == true
#          select1   = MAF .< vec(p) .< 1-MAF
#          select2   = vec(var(genotypes,dims=1)) .!= 0
#          select    = select1 .& select2
#          genotypes = genotypes[:,select]
#          p         = p[:,select]
#          markerID  = markerID[select]
#          printstyled("$(sum(1 .- select)) loci which are fixed or have minor allele frequency < $MAF are removed.\n",bold=true)
#     end
#     nObs,nMarkers = size(genotypes)       #number of individuals and molecular markers
#     sum2pq        = (2*p*(1 .- p)')[1,1]  #∑2pq

#     # #a kernel / relationship matrix or genotypes are provided in GBLUP
#     # if method == "GBLUP"
#     #     if !isGRM #calculate the relationship matrix from the genotype covariate matrix
#     #         genotypes  = genotypes ./ sqrt.(2*p.*(1 .- p))
#     #         genotypes  = (genotypes*genotypes'+ I*data_type(0.00001))/nMarkers
#     #         println("A genomic relationship matrix is computed from genotypes.")
#     #         isGRM  = true
#     #     end
#     #     add_small_value_count = 0
#     #     while isposdef(genotypes) == false
#     #         println("The relationship matrix is not positive definite. A very small number is added to the diagonal.")
#     #         genotypes = genotypes + I*data_type(0.00001)
#     #         add_small_value_count += 1
#     #         if add_small_value_count > 10
#     #             error("Please provide a positive-definite relationship matrix.")
#     #         end
#     #     end
#     # end

#     genotypes     = nnmmGenotypes(obsID,markerID,nObs,nMarkers,p,sum2pq,center,genotypes,isGRM)

#     # genotypes.G                = Variance(G_is_marker_variance ? G : false, df,false,estimate_variance,estimate_scale,constraint) 
#     # genotypes.genetic_variance = Variance(G_is_marker_variance ? false : G, df,false,estimate_variance,estimate_scale,constraint) 
    
#     # genotypes.method     = method
#     # genotypes.estimatePi = estimatePi
#     # genotypes.π          = Pi

#     writedlm("IDs_for_individuals_with_genotypes.txt",genotypes.obsID)
#     println("Genotype information:")
#     println("#markers: ",(isGRM ? 0 : size(genotypes.data,2)),"; #individuals: ",size(genotypes.data,1))

#     #starting values for marker effects
#     if starting_value != false
#         printstyled("Starting values are provided for marker efffects. The order of starting values should be\n",
#         "markers for each traits (all markers for trait 1 then all markers for trait 2...)\n",bold=false,color=:green)
#         nsol = (method != "GBLUP" ? nMarkers : nObs)
#         if length(starting_value)%nsol != 0 #work for both single and multiple traits
#             error("length of starting values is wrong.")
#         end
#         genotypes.α = starting_value
#     end
#     return genotypes
# end

function read_omics(file::AbstractString, G = false;
                       ## method:
                       method = "BayesC", Pi = 0.0, estimatePi = true, 
                       ## variance:
                       G_is_marker_variance = false, df = 4.0,
                       estimate_variance = true, estimate_scale = false,
                       constraint = false, #for multi-trait only, constraint=true means no genetic covariance among traits
                       # format:
                       separator = ',', header = true, double_precision = false,
                       ## quality control:
                       quality_control = true, missing_value = "NA",
                       ## others:
                       center = false, starting_value = false)
    #Define data precision
    data_type = double_precision ? Float64 : Float32 #works for both genotypes and GRM
    #Read the genotype file
    #print info about the file format
    printstyled("The delimiter in ", split(file, ['/', '\\'])[end], " is \'", separator, "\'. ", bold=false, color=:green)
    printstyled("The header (omics IDs) is ", (header ? "provided" : "not provided"), " in ", split(file, ['/', '\\'])[end], ".\n", bold=false, color=:green)
    data = CSV.read(file,DataFrame,delim = separator, header=header,missingstrings=missing_value)

    obsID     = map(string,data[!,1]) # individual ID
    omicsID   = names(data)[2:end]    #1st column is individual ID
    omics      = Matrix(data[!,2:end]) #memory usage is doubled here
    
    #clean memory
    data = nothing
    GC.gc()

    #preliminary summary of genotype
    nObs,nOmics = size(omics)       #number of individuals and molecular markers

    #Naive Quality Control 1, warn omics with large values
    if quality_control == true
        for omicsi in eachcol(omics)
            missing_obs        = findall(x->x==float(missing_value),omicsi)  ## How CSV.read deal with `missingstrings`?? are they already `missing`?
            omicsi[missing_obs] .= mean(omicsi[nonmissing_obs])
            std_omicsi=std(omicsi)
            if findfirst(x->(x>3*std_omicsi || x<3*std_omicsi),omicsi) != nothing 
                @warn "omics out of the 3*std are found."
            end
        end

    end

    omicsMeans   = center==true ? center!(omics) : mean(omics,dims=1) #centering omics or not

    #Naive Quality Control 2, fixed omics
    if quality_control == true
         select   = vec(var(omics,dims=1)) .!= 0
         omics = omics[:,select]
        #  p         = p[:,select]
         omicsID  = omicsID[select]
         printstyled("$(sum(1 .- select)) omics which are fixed are removed.\n",bold=true)
    end
    nObs,nOmics = size(omics)       #number of individuals and omics
 
    omics     = Omics(obsID,omicsID,nObs,nOmics,center,omics)

    writedlm("IDs_for_individuals_with_omics.txt",omics.obsID)
    println("Omics information:")
    println("#omics: ", omics.nFeatures,"; #individuals: ", omics.nObs) #use omics.nFeatures instead of nOmics

    #starting values for omics effects
    if starting_value != false
        printstyled("Starting values are provided for omics effects. The order of starting values should be\n",
        "omics for each traits (all omics for trait 1 then all omics for trait 2...)\n",bold=false,color=:green)
        nsol = nMarkers
        if length(starting_value)%nsol != 0 #work for both single and multiple traits
            error("length of starting values is wrong.")
        end
        omics.α = starting_value
    end
    return omics
end

function read_phenotypes(file::AbstractString;
                       ## format:
                       separator = ',', header = true, double_precision = false, missing_value = "NA")
    #Define data precision
    data_type = double_precision ? Float64 : Float32 #works for both genotypes and GRM
    #print info about the file format
    printstyled("The delimiter in ", split(file, ['/', '\\'])[end], " is \'", separator, "\'. ", bold=false, color=:green)
    printstyled("The header (pheno IDs) is ", (header ? "provided" : "not provided"), " in ", split(file, ['/', '\\'])[end], ".\n", bold=false, color=:green)
    #read a large pheno file
    data      = CSV.read(file,DataFrame,delim = separator,header=header,missingstrings=missing_value)
    data[!,1]  = convert.(String, data[!,1]) #convert the first column to string
    obsID     = data[!,1]
    phenoID    = names(data)[2:end]
    pheno      = Matrix(data[!,2:end])


    #preliminary summary of genotype
    nObs,nPheno = size(pheno)       #number of individuals and molecular markers

    pheno  = Phenotypes(obsID,phenoID,nObs,nPheno,pheno)

    writedlm("IDs_for_individuals_with_phenotypes.txt",pheno.obsID)
    println("Pheno information:")
    println("#pheno: ", nPheno,"; #individuals: ", pheno.nObs)

    return pheno
end

