include("GWAS_deprecated.jl")

"""
    GWAS(marker_effects_file;header=true)

Compute the model frequency for each marker (the probability the marker is included in the model) using samples of marker effects stored in **marker_effects_file**.
"""
function GWAS(marker_effects_file::AbstractString...;header=true)
    is_nnbayes = length(marker_effects_file)==1 ? false : true
    if is_nnbayes && header==false
        error("Header is required for all files.")
    end

    println("Compute the model frequency for each marker (the probability the marker is included in the model).")
    file = marker_effects_file

    #read all files
    markerID=[]
    samples=[]
    for i in 1:length(marker_effects_file)
        if header==true
            sample_i,markerID_i =readdlm(file[i],',',header=true)
        else
            sample_i = readdlm(file[i],',')
            markerID_i = 1:size(sample_i,2)
        end
        append!(samples,[sample_i])
        append!(markerID,[markerID_i])
    end

    if is_nnbayes #NNBayes
        num_hidden_nodes = length(marker_effects_file)
        println("NN-Bayes:")
        println(" - MCMC samples for marker effects on $num_hidden_nodes hidden nodes were loaded.")
        if all(id->id==markerID[1], markerID)  #id in all files are the same -> fully-connected
            println(" - fully-connected neural networks.")
            markerID = markerID[1]                    # only keep id from one file
            samples  = map(x -> x.!=0, samples)  # change non-zero effects to 1
            samples  = sum(samples)              # indicate elements that are non-zero on at least one of the hidden nodes
        else   #different id in different files -> partial-connected
            println(" - partial-connected neural networks.")
            markerID = hcat(markerID...)  #merge id together
            samples  = hcat(samples...)   #merge samples together
        end
    else #conventional linear models
        samples = samples[1]
        markerID = markerID[1]
    end

    modelfrequency = vec(mean(samples .!= 0.0,dims=1))
    out = DataFrame(marker_ID = vec(markerID), modelfrequency = modelfrequency)
    return out
end

"""
    GWAS(model,map_file,marker_effects_file...;
         window_size = "1 Mb",sliding_window = false,
         GWAS = true, threshold = 0.001,
         genetic_correlation = false,
         header = true)

run genomic window-based GWAS

* MCMC samples of marker effects are stored in **marker_effects_file** with delimiter ','.
* **model** is either the model::MME used in analysis or the genotype cavariate matrix M::Array
* **map_file** has the (sorted) marker position information with delimiter ','. If the map file is not provided,
  i.e., **map_file**=`false`, a fake map file will be generated with **window_size** markers in each 1 Mb window, and
  each 1 Mb window will be tested.
* If two **marker_effects_file** are provided, and **genetic_correlation** = true, genomic correlation for each window is calculated.
* Statistics are computed for nonoverlapping windows of size *window_size* by default.
  If **sliding_window** = true, those for overlapping sliding windows are calculated.
* map file format:

```
markerID,chromosome,position
m1,1,16977
m2,1,434311
m3,1,1025513
m4,2,70350
m5,2,101135
```

"""
function GWAS(mme,map_file,marker_effects_file::AbstractString...;
              #NNBayes
              weights_NN_file = false,  #MCMC samples for NN weights, e.g., "MCMC_samples_neural_networks_bias_and_weights.txt"
              ebv_file = false,         #MCMC samples for ebv of each individual, e.g., "MCMC_samples_EBV_NonLinear.txt"
              nonlinear_function = false,
              #window
              window_size = "1 Mb",sliding_window = false,
              #GWAS
              GWAS = true, threshold = 0.001,
              #genetic correlation
              genetic_correlation = false,
              #local EBV
              local_EBV = false,
              #misc
              header = true, output_winVarProps = false)
    #NN-Bayes
    is_nnbayes = weights_NN_file!=false && ebv_file!=false && nonlinear_function!=false ? true : false
    if is_nnbayes
        #covert nonlinear function
        nonlinear_function = nnbayes_activation(nonlinear_function)

        #read data
        marker_effects=[]
        markerID=[]
        num_hidden_nodes = length(marker_effects_file)
        for i in 1:num_hidden_nodes
            marker_effects_i,markerID_i = readdlm(marker_effects_file[i],',',header=true)
            append!(marker_effects,[marker_effects_i])
            append!(markerID,[markerID_i])
        end
        weights_NN = readdlm(weights_NN_file,',')
        ebv = readdlm(ebv_file,header=true,',')[1]

        #print info
        println("NN-Bayes:")
        println(" - MCMC samples for marker effects on $num_hidden_nodes hidden nodes were loaded.")
        if all(id->id==markerID[1], markerID)  #id in all files are the same -> fully-connected
            println(" - fully-connected neural networks.")
        else   #different id in different files -> partial-connected
            error("WPPA for partial-connected neural networks is not supported!")
        end
    end

    #get genotype
    X = (typeof(mme) <: Array ? mme : mme.M[1].output_genotypes)
    nmarkers=size(X,2)

    if typeof(window_size) == String
        if split(window_size)[2] != "Mb"
            error("The format for window_size is \"1 Mb\".")
        end
    end

    if map_file == false && typeof(window_size) <: Integer
        println("The map file is not provided. A fake map file is generated with $window_size markers in each 1 Mb window.")
<<<<<<< HEAD
=======
        nmarkers=length(readdlm(marker_effects_file[1],',',header=true)[2])
>>>>>>> master
        mapfile = DataFrame(markerID  =1:nmarkers,
                            chromosome=fill(1,nmarkers),
                            position  =Int.(floor.(1:(1_000_000/window_size):nmarkers*(1_000_000/window_size))))
        CSV.write("mapfile.temp",mapfile)
        map_file, window_size = "mapfile.temp", "1 Mb"
    end

    window_size_bp = map(Int64,parse(Float64,split(window_size)[1])*1_000_000)
    mapfile = (header == true ? readdlm(map_file,',',header=true)[1] : readdlm(map_file,','))
    chr     = map(string,mapfile[:,2])
    pos     = map(Int64,mapfile[:,3])

    window_size_nSNPs   = Array{Int64,1}()  #save number of markers in ith window for all windows
    window_chr          = Array{String,1}() #1
    window_pos_start    = Array{Int64,1}()  #1_000_000
    window_pos_end      = Array{Int64,1}()  #2_000_000
    window_snp_start    = Array{Int64,1}()  #1_314_314
    window_snp_end      = Array{Int64,1}()  #1_999_003
    window_column_start = Array{Int64,1}()  #101
    window_column_end   = Array{Int64,1}()  #200

    index_start = 1
    for i in unique(chr)
      pos_on_chri     = pos[chr.== i] #assume chr and pos are sorted
      if sliding_window == false
          nwindow_on_chri = ceil(Int64,pos_on_chri[end]/window_size_bp)
      else
          nwindow_on_chri = findfirst(x -> x >= pos_on_chri[end] - window_size_bp, pos_on_chri)
      end

      for j in 1: nwindow_on_chri
        if sliding_window == false
            thisstart = window_size_bp*(j-1)
        else
            thisstart = pos_on_chri[j]
        end
        thisend  = thisstart + window_size_bp
        snps_window = thisstart .<= pos_on_chri .< thisend
        snps_window_sizej = sum(snps_window)
        #empty windows exist in non-sliding window; no empty window in sliding windows
        #empty windows were deleted
        if snps_window_sizej!=0
            push!(window_snp_start,pos_on_chri[findfirst(snps_window)])
            push!(window_snp_end,pos_on_chri[findlast(snps_window)])
            push!(window_column_start,index_start)
            push!(window_column_end,index_start+snps_window_sizej-1)
            push!(window_chr,i)
            push!(window_pos_start,thisstart)
            push!(window_pos_end,thisend)
            push!(window_size_nSNPs,snps_window_sizej)
        end
        if sliding_window == false
            index_start += snps_window_sizej
        else
            index_start += 1
        end
      end
    end

    out=[]
    if GWAS == true
        println("Compute the posterior probability of association of the genomic window that explains more than ",threshold," of the total genetic variance.")
        niter = is_nnbayes ? 1 : length(marker_effects_file)
        for i in 1:niter
            #using marker effect files
            output            = is_nnbayes ?  marker_effects : readdlm(marker_effects_file[i],',',header=true)[1]
            nsamples          = is_nnbayes ?  size(output[1],1) : size(output,1)
            nWindows          = length(window_size_nSNPs)
            winVarProps       = zeros(nsamples,nWindows)
            winVar            = zeros(nsamples,nWindows)
            #window_mrk_start ID and window_mrk_end ID are not provided now
<<<<<<< HEAD

=======
            X = (typeof(mme) <: Array ? mme : mme.M[1].output_genotypes)
            if local_EBV==true
                nind     = size(X,1)
                localEBV = zeros(nind,nWindows)
            end
>>>>>>> master
            @showprogress "running GWAS..." for i=1:nsamples
                if is_nnbayes
                    α = hcat(map(x -> x[i,:], output)...)
                    w = weights_NN[i,:]
                    ebv_i = ebv[i,:]
                    genVar = var(ebv_i)
                else
                    α = output[i,:]
                    genVar = var(X*α)
                end

                for winj = 1:length(window_column_start)
                  wStart = window_column_start[winj]
                  wEnd   = window_column_end[winj]
                  if is_nnbayes
                      BV_winj = nonlinear_function.(X[:,wStart:wEnd]*α[wStart:wEnd,:])*w[2:end]
                  else #conventional
                      BV_winj = X[:,wStart:wEnd]*α[wStart:wEnd]
                  end
                  var_winj = var(BV_winj)
                  winVar[i,winj]      = var_winj
                  winVarProps[i,winj] = var_winj/genVar
                  if local_EBV == true
                      localEBV[:,winj] += (BV_winj - localEBV[:,winj])/i
                  end
                end
            end
            if local_EBV == true
                df= DataFrame(localEBV)
                rename!(df,Symbol.("w".*string.(1:nWindows)))
                CSV.write("localEBV"*string(i)*".txt",hcat(DataFrame(ID=mme.output_ID),df))
            end
            winVarProps[isnan.(winVarProps)] .= 0.0 #replace NaN caused by situations no markers are included in the model
            WPPA, prop_genvar = vec(mean(winVarProps .> threshold,dims=1)), vec(mean(winVarProps,dims=1))
            prop_genvar = round.(prop_genvar*100,digits=2)
            winVarmean = vec(mean(winVar,dims=1))
            winVarstd  = vec(std(winVar,dims=1))

            srtIndx = sortperm(WPPA,rev=true)

            outi = DataFrame(trait  = fill(i,length(WPPA))[srtIndx],
                            window = (1:length(WPPA))[srtIndx],
                            chr    = window_chr[srtIndx],
                            wStart = window_pos_start[srtIndx],
                            wEnd   = window_pos_end[srtIndx],
                            start_SNP = window_snp_start[srtIndx],
                            end_SNP   = window_snp_end[srtIndx],
                            numSNP  = window_size_nSNPs[srtIndx],
                            estimateGenVar  = winVarmean[srtIndx],
                            stdGenVar     = winVarstd[srtIndx],
                            prGenVar = prop_genvar[srtIndx],
                            WPPA     = WPPA[srtIndx],
                            PPA_t  = cumsum(WPPA[srtIndx]) ./ (1:length(WPPA)))
             push!(out,outi)
             outfile="GWAS_"*replace(string.(marker_effects_file[i]),"/"=>"_")
             CSV.write(outfile, outi)
        end
    end
    if genetic_correlation == true && length(marker_effects_file) ==2
            #using marker effect files
            output1           = readdlm(marker_effects_file[1],',',header=true)[1]
            output2           = readdlm(marker_effects_file[2],',',header=true)[1]
            nsamples,nmarkers = size(output1)
            nWindows          = length(window_size_nSNPs)
            gcov              = zeros(nsamples,nWindows)
            gcor              = zeros(nsamples,nWindows)
            #window_mrk_start ID and window_mrk_end ID are not provided now
            X = (typeof(mme) <: Array ? mme : mme.M[1].output_genotypes)
            @showprogress "calculating genomic correlation..." for i=1:nsamples
                α1 = output1[i,:]
                α2 = output2[i,:]
                for winj = 1:length(window_column_start)
                  wStart = window_column_start[winj]
                  wEnd   = window_column_end[winj]
                  Xwinj  = X[:,wStart:wEnd]
                  BV_winj1= Xwinj*α1[wStart:wEnd]
                  BV_winj2= Xwinj*α2[wStart:wEnd]
                  gcov[i,winj] = cov(BV_winj1,BV_winj2)
                  gcor[i,winj] = cor(BV_winj1,BV_winj2)
                end
            end
            gcov[isnan.(gcov)] .= 0.0
            gcor[isnan.(gcor)] .= 0.0
            gcovmean,gcovstd = vec(mean(gcov,dims=1)),vec(std(gcov,dims=1))
            gcormean,gcorstd = vec(mean(gcor,dims=1)),vec(std(gcor,dims=1))
            outi = DataFrame(trait  = fill("cor(t1,t2)",length(gcormean)),
                            window = (1:length(gcormean)),
                            chr    = window_chr,
                            wStart = window_pos_start,
                            wEnd   = window_pos_end,
                            start_SNP = window_snp_start,
                            end_SNP   = window_snp_end,
                            numSNP   = window_size_nSNPs,
                            estimate_cov = gcovmean,
                            std_cov      = gcovstd,
                            estimate_cor = gcormean,
                            std_cor      = gcorstd)
             push!(out,outi)
             outfile="GWAS_"*replace(string.(marker_effects_file[1]),"/"=>"_")*"_"*replace(string.(marker_effects_file[2]),"/"=>"_")
             CSV.write(outfile, outi)
    end
    return output_winVarProps ? (Tuple(out),winVarProps) : Tuple(out)
end

function GWAS(marker_effects_file::AbstractString,map_file::AbstractString,mme;
             header=true,
             window_size="1 Mb",
             threshold=0.001,
             sliding_window = false,
             output_winVarProps=false,
             local_EBV = false)
     GWAS(mme,map_file,marker_effects_file,
                  header=header,
                  window_size=window_size,
                  threshold=threshold,
                  sliding_window = sliding_window,
                  output_winVarProps=output_winVarProps,
                  local_EBV = local_EBV)
end


export GWAS


# below function is to define the activation function for neural network
function nnbayes_activation(activation_function)
      if activation_function == "tanh"
         mytanh(x) = tanh(x)
         return mytanh
      elseif activation_function == "sigmoid"
         mysigmoid(x) = 1/(1+exp(-x))
         return mysigmoid
      elseif activation_function == "relu"
         myrelu(x) = max(0, x)
         return myrelu
      elseif activation_function == "leakyrelu"
         myleakyrelu(x) = max(0.01x, x)
         return myleakyrelu
      elseif activation_function == "linear"
         mylinear(x) = x
         return mylinear
      else
          error("invalid actication function")
      end
end
