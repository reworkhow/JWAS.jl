"""
    get_column_ref(X::Vector{T})

    To obtain a vector of views (alias/pointer) for each column of the input matrix WITHOUT COPYING the underlying data. 
    input:  a matrix X
    output: a vector containing views of each column of the input matrix X
"""
function get_column_ref(X::AbstractArray{T,2}) where T <: Any
      return [view(X, :, i) for i in 1:size(X, 2)] # Create a vector of views for each column of the matrix X
end

"""
    This function centers columns of the input matrix X by subtracting their means along each column. The function operates in-place by modifying the original matrix X.

    Input:
    - X::AbstractMatrix: a matrix to be centered

    Output:
    - col_means::Vector: a vector of mean values for each column in the original matrix, computed before centering.
"""
function center!(X::AbstractArray{T,2}) where T <: Any
    col_means = mean(X, dims=1) # Calculate the means in each column
    X .-= col_means  # Subtract column means from X, 
                     # Note: in-place, using broadcasting, i.e. X[i,j] -= col_means[j]
    return col_means # Return a row vector of column means (size 1 by ncol)
end

function getXpRinvX(X, Rinv)
    XpRinvX = [((X[:,i].*Rinv)'X[:,i]) for i=1:size(X,2)]
    return XpRinvX
end

function getXpRinvX(X)
    XpRinvX = [dot(X[:,i],X[:,i]) for i=1:size(X,2)]
    return XpRinvX
end

"""
    is_unit_weights(Rinv)

Return `true` if all residual weights are exactly one.
"""
function is_unit_weights(Rinv::AbstractVector)
    oneval = one(eltype(Rinv))
    @inbounds for i in eachindex(Rinv)
        if Rinv[i] != oneval
            return false
        end
    end
    return true
end

"""
    block_rhs!(rhs, Xblock, y, Rinv, unit_weights)

Compute block RHS `Xblock' * Diagonal(Rinv) * y` in-place without storing
`Xblock' * Diagonal(Rinv)` as a persistent matrix.
"""
function block_rhs!(rhs::AbstractVector,
                    Xblock::AbstractMatrix,
                    y::AbstractVector,
                    Rinv::AbstractVector,
                    unit_weights::Bool)
    if unit_weights
        mul!(rhs,transpose(Xblock),y)
        return rhs
    end

    @inbounds for j in 1:size(Xblock,2)
        xj = view(Xblock,:,j)
        acc = zero(eltype(rhs))
        @simd for i in eachindex(y)
            acc += xj[i]*Rinv[i]*y[i]
        end
        rhs[j] = acc
    end
    return rhs
end

function get_column_blocks_ref(X,fast_blocks=false)
    xArray = Array{typeof(X)}(undef,length(fast_blocks))
    @showprogress "building prerequisite matrices ..." for i in 1:length(fast_blocks)
        pos_start = fast_blocks[i]
        pos_end   = (i != length(fast_blocks)) ? (fast_blocks[i+1]-1) : size(X,2)
        xArray[i] = view(X,:,pos_start:pos_end)
    end
    return xArray
end

"""
    estimate_marker_memory(nObs, nMarkers;
                           element_bytes,
                           has_nonunit_weights=false,
                           block_starts=false)

Estimate major marker-path memory components in bytes.
"""
function estimate_marker_memory(nObs::Integer,nMarkers::Integer;
                                element_bytes::Integer,
                                has_nonunit_weights::Bool=false,
                                block_starts=false)
    if nObs < 0 || nMarkers < 0
        error("nObs and nMarkers must be non-negative.")
    end
    if element_bytes <= 0
        error("element_bytes must be a positive integer.")
    end

    # Use wider integers so large N*P products do not overflow on 64-bit Int.
    nObs128        = Int128(nObs)
    nMarkers128    = Int128(nMarkers)
    bytes_per_val  = Int128(element_bytes)

    nvals_x        = nObs128*nMarkers128
    bytes_X        = nvals_x*bytes_per_val
    bytes_xpRinvx  = nMarkers128*bytes_per_val
    bytes_xRinv    = has_nonunit_weights ? bytes_X : Int128(0)

    is_block_mode  = block_starts != false
    # XRinvArray is no longer persisted in block mode; block RHS is computed on demand.
    bytes_XRinv    = Int128(0)
    block_sq_sum   = Int128(0)
    if is_block_mode
        starts = collect(Int,block_starts)
        nblocks = length(starts)
        for i in 1:nblocks
            block_start = starts[i]
            block_end   = (i == nblocks) ? nMarkers : (starts[i+1]-1)
            block_size  = max(0,block_end-block_start+1)
            # XpRinvX stores per-block Gram matrices, so memory scales with sum(s_i^2).
            block_sq_sum += Int128(block_size)*Int128(block_size)
        end
    end
    bytes_XpRinvX  = block_sq_sum*bytes_per_val
    bytes_total    = bytes_X + bytes_xRinv + bytes_XRinv + bytes_XpRinvX + bytes_xpRinvx

    return (
        bytes_X = bytes_X,
        bytes_xRinvArray = bytes_xRinv,
        bytes_XRinvArray = bytes_XRinv,
        bytes_XpRinvX = bytes_XpRinvX,
        bytes_xpRinvx = bytes_xpRinvx,
        bytes_total = bytes_total
    )
end

"""
    format_bytes_human(bytes)

Format byte counts into a compact human-readable string.
"""
function format_bytes_human(bytes::Integer)
    if bytes < 0
        error("bytes must be non-negative.")
    end
    units = ("B","KiB","MiB","GiB","TiB","PiB","EiB")
    value = Float64(bytes)
    unit_index = 1
    while value >= 1024 && unit_index < length(units)
        value /= 1024
        unit_index += 1
    end
    return @sprintf("%.2f %s",value,units[unit_index])
end

"""
    check_marker_memory_guard!(; mode, ratio, estimated_bytes, total_memory_bytes, context_string)

Apply guard policy for estimated marker memory usage.
"""
function check_marker_memory_guard!(; mode::Symbol=:error,
                                     ratio::Real=0.8,
                                     estimated_bytes::Integer,
                                     total_memory_bytes::Integer,
                                     context_string::AbstractString="")
    if !(mode in (:error,:warn,:off))
        error("memory_guard must be one of :error, :warn, or :off.")
    end
    if !(0 < ratio <= 1)
        error("memory_guard_ratio must satisfy 0 < memory_guard_ratio <= 1.")
    end
    if estimated_bytes < 0 || total_memory_bytes <= 0
        error("estimated_bytes must be >= 0 and total_memory_bytes must be > 0.")
    end

    if mode == :off
        return :skipped
    end

    # Guard compares estimated marker memory against a configurable RAM fraction.
    threshold = Int128(floor(Float64(total_memory_bytes)*Float64(ratio)))
    estimated = Int128(estimated_bytes)
    if estimated <= threshold
        return :ok
    end

    msg = "Estimated marker memory usage exceeds configured guard threshold.\n"*
          "context: $(context_string)\n"*
          "estimated: $(format_bytes_human(estimated))\n"*
          "threshold ($(ratio*100)% of RAM): $(format_bytes_human(threshold))\n"*
          "system RAM: $(format_bytes_human(total_memory_bytes))\n"*
          "Set memory_guard=:warn or :off to override, or reduce model/data size."

    if mode == :warn
        @warn msg
        return :warned
    end
    error(msg)
end

mutable struct GibbsMats
    X::Union{Array{Float64,2},Array{Float32,2}}
    nrows::Int64
    ncols::Int64
    xArray     #do not declare type because it is a vector of views
    xRinvArray #do not declare type because it is a vector of views
    xpRinvx::Union{Array{Float64,1},Array{Float32,1}}
    XArray
    XRinvArray
    XpRinvX
    function GibbsMats(X::Union{Array{Float64,2},Array{Float32,2}},Rinv;fast_blocks=false)
        #single-locus, adjusting y, approach
        nrows,ncols = size(X)
        xArray      = get_column_ref(X)
        unit_weights = is_unit_weights(Rinv)
        xpRinvx     = unit_weights ? getXpRinvX(X) : getXpRinvX(X,Rinv)
        if unit_weights
            xRinvArray = xArray  #avoid using extra memory for xRinvArray
        else
            xRinvArray = [x.*Rinv for x in xArray]
        end
        #block, adjusting rhs and/or y, approach
        if fast_blocks != false
            XArray = get_column_blocks_ref(X,fast_blocks)
            XRinvArray = false # Save memory by not persisting X'Rinv; block RHS is computed on demand.
            if unit_weights
                XpRinvX = @showprogress "building prerequisite matrices ..." [Xblock' * Xblock for Xblock in XArray]
            else
                Rdiag = Diagonal(Rinv)
                XpRinvX = @showprogress "building prerequisite matrices ..." [Xblock' * Rdiag * Xblock for Xblock in XArray]
            end
        else
            XArray = XRinvArray = XpRinvX = false
        end
        new(X,nrows,ncols,
            xArray,xRinvArray,xpRinvx,
            XArray,XRinvArray,XpRinvX)
    end
end

################################################################################
#align genotypes with phenotyps given IDs #MAY USE TOO MUCH MEMORY
################################################################################
#input: mme.M.genotypes (1) complete genomic data  -> all genotyped individuals
#                       (2)* incomplete genomic data -> all pedigree individuals
#output:mme.M.genotypes (1) subset
#
# mme.M.genotypes is modified to same IDs in phenotypes for coding convenience
#
#Note that Aligning genotypes with phenotypes and output_ID in single-step analysis
#(incomplete genomic data) has been moved to SSBR.jl file.
function align_genotypes(mme::MME,output_heritability=false,single_step_analysis=false)
    if single_step_analysis==false
        if mme.output_ID != 0
            for Mi in mme.M
                if mme.output_ID != Mi.obsID
                    Zo  = map(Float32,mkmat_incidence_factor(mme.output_ID,Mi.obsID))
                    Mi.output_genotypes =  Zo*Mi.genotypes
                else
                    Mi.output_genotypes =  Mi.genotypes #reference, not copy
                end
                if Mi.isGRM #relationship matrix is provided
                    Z  = map(Float32,mkmat_incidence_factor(mme.obsID,Mi.obsID))
                    Mi.output_genotypes = Mi.output_genotypes*Z'
                end
            end
        end
        #***************************************************************************
        #Align genotypes with phenotypes
        #Change mme.M.genotypes to phenotyped individuals to accommodate
        #individuals with repeated records or individuals without records
        #
        #**********CENTERING?*******************************************************
        if mme.obsID != mme.M[1].obsID
            for Mi in mme.M
                Z  = map(Float32,mkmat_incidence_factor(mme.obsID,Mi.obsID))
                genotypes = Z*Mi.genotypes
                if Mi.isGRM #relationship matrix is provided
                    genotypes = genotypes*Z'
                end
                Mi.genotypes = genotypes
                Mi.obsID     = mme.obsID
                Mi.nObs      = length(mme.obsID)
            end
        end
    end
end

"""
mkmat_incidence_factor(yID::Vector, uID::Vector)
    create an incidence matrix Z to reorder uID to yID by yID = Z*uID.
    input: 
        - yID: a vector containing the desired order of IDs
        - uID: a vector containing the original order of IDs
    output: 
        - Z: a sparse matrix representing the incidence relationship between yID and uID (yID = Z*uID)
"""
function mkmat_incidence_factor(yID::Vector, uID::Vector)
    Z = spzeros(length(yID), length(uID)) # initialize a sparse matrix Z with the dimensions of yID and uID
    uIDdict = Dict(id => index for (index, id) in enumerate(uID))  # create a uID dictionary to map uID elements to their indices
    
    # iterate over yID elements to populate Z matrix
    for (rowi, id) in enumerate(yID)
        if haskey(uIDdict, id) # check if the current yID element exists in uID dictionary
            index = uIDdict[id]
        else
            error("$id is not found!") # error if the id is not found in the uID dictionary
        end
        Z[rowi, index] = 1 # set the corresponding element in the Z matrix to 1
    end
    return Z
end

################################################################################
#Set Marker Effcets Hyperparameters: Variances and Pi
################################################################################
function set_marker_hyperparameters_variances_and_pi(mme::MME)
    if mme.M != 0
        for Mi in mme.M
            #(1) Pi in multi-trait analysis
            if (mme.nModels !=1 || mme.MCMCinfo.RRM != false) && Mi.π==0.0 && !(Mi.method in ["RR-BLUP","BayesL"])
                println()
                printstyled("Pi (Π) is not provided.\n",bold=false)
                printstyled("Pi (Π) is generated assuming all markers have effects on all traits.\n",bold=false)
                mykey=Array{Float64}(undef,0)
                ntraits=mme.nModels
                if mme.MCMCinfo.RRM != false
                    ntraits = size(mme.MCMCinfo.RRM,2)
                end
                Pi=Dict{Array{Float64,1},Float64}()
                #for i in [ bin(n,ntraits) for n in 0:2^ntraits-1 ] `bin(n, pad)` is deprecated, use `string(n, base=2, pad=pad)
                for i in [ string(n,base=2,pad=ntraits) for n in 0:2^ntraits-1 ]
                    Pi[parse.(Float64, split(i,""))]=0.0
                end
                Pi[ones(ntraits)]=1.0
                Mi.π = Pi
            end
            #(2) marker effect variances
            if Mi.G.val == false
                if Mi.method!="GBLUP"
                    genetic2marker(Mi,Mi.π)
                    println()
                    if mme.nModels != 1 || mme.MCMCinfo.RRM != false
                      if !isposdef(Mi.G.val) #also work for scalar
                        error("Marker effects covariance matrix is not postive definite! Please modify the argument: Pi.")
                      end
                      println("The prior for marker effects covariance matrix is calculated from genetic covariance matrix and Π.")
                      println("The mean of the prior for the marker effects covariance matrix is:")
                      if mme.MCMCinfo.printout_model_info==true
                          Base.print_matrix(stdout,round.(Mi.G.val,digits=6))
                      end
                    else
                      if !isposdef(Mi.G.val) #positive scalar (>0)
                        error("Marker effects variance is negative!")
                      end
                      println("The prior for marker effects variance is calculated from the genetic variance and π.")
                      print("The mean of the prior for the marker effects variance is: ")
                      print(round.(Mi.G.val,digits=6))
                    end
                    print("\n\n\n")
                elseif Mi.method == "GBLUP"
                    Mi.G.val  = Mi.genetic_variance.val
                end
            end
            #(3) scale parameter for marker effect variance
            if Mi.ntraits == 1 && mme.MCMCinfo.RRM == false
                Mi.G.scale = Mi.G.val*(Mi.G.df-Mi.ntraits-1)/Mi.G.df
            else
                Mi.G.scale = Mi.G.val*(Mi.G.df-Mi.ntraits-1)
            end
        end
    end
end

function genetic2marker(M::Genotypes,Pi::Dict)
  ntraits = size(M.genetic_variance.val,1)
  denom   = zeros(ntraits,ntraits)
  for i in 1:ntraits
    for j in i:ntraits
      pi_selected = filter(d->d.first[i]==1.0 && d.first[j]==1.0,Pi)

      denom[i,j] = M.sum2pq*sum(values(pi_selected))
      denom[j,i] = denom[i,j]
    end
  end
  M.G.val = M.genetic_variance.val ./ denom
end

function genetic2marker(M::Genotypes,π::Float64)
    M.G.val = M.genetic_variance.val/((1-π)*M.sum2pq)
end
