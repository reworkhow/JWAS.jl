include("GWAS_deprecated.jl")

"""
    GWAS(marker_effects_file;header=true)

Compute the model frequency for each marker (the probability the marker is included in the model) using samples of marker effects stored in **marker_effects_file**.
"""
function GWAS(marker_effects_file;header=true)
    println("Compute the model frequency for each marker (the probability the marker is included in the model).")
    file = marker_effects_file
    if header==true
        samples,markerID =readdlm(file,',',header=true)
    else
        samples=readdlm(file,',')
        markerID = 1:size(samples,2)
    end
    modelfrequency = vec(mean(samples .!= 0.0,dims=1))
    out = DataFrame(marker_ID = vec(markerID), modelfrequency = modelfrequency)
    return out
end

#Chen, C., Steibel, J. P., & Tempelman, R. J. (2017). Genome-Wide Association
#Analyses Based on Broadly Different Specifications for Prior Distributions,
#Genomic Windows, and Estimation Methods. Genetics, 206(4), 1791–1806.
#Also, as per Moser et al. (2015), two different within-chromosome starting positions
#(starting at location 0 or 0.25 Mb for window sizes 0.5, starting at 0 or 0.5 Mb
#location for window sizes 1 Mb, starting at 0 or 1 Mb location for window sizes 2 Mb,
#and starting at 0 or 1.5 Mb location for window sizes 3 Mb) for each chromosome were
#chosen to partly counteract the chance effect of different LD patterns being associated
#with nonoverlapping windows. Finally, adaptive window sizes based on clustering SNP
#by LD r2 were also determined using the BALD R package (Dehman and Neuvial 2015),
#using the procedure described by Dehman et al. (2015).
"""
    GWAS(model,map_file,marker_effects_file...;
         window_size = "1 Mb",sliding_window = false,
         GWAS = true, threshold = 0.001,
         genetic_correlation = false,
         header = true)

run genomic window-based GWAS

* MCMC samples of marker effects are stored in **marker_effects_file** with delimiter ','.
* **model** is either the model::MME used in analysis or the genotype cavariate matrix M::Array
* **map_file** has the (sorted) marker position information with delimiter ','. If the map file is not provided, a fake map file
  will be generated with 100 markers in each 1 Mb window.
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
function GWAS(mme,map_file::AbstractString,marker_effects_file::AbstractString...;
              #window
              window_size = "1 Mb",sliding_window = false,
              #GWAS
              GWAS = true, threshold = 0.001,
              #genetic correlation
              genetic_correlation = false,
              #misc
              header = true, output_winVarProps = false)

    if split(window_size)[2] != "Mb"
        error("The format for window_size is \"1 Mb\".")
    end
    if map_file == false
        println("The map file is not provided. A fake map file is generated with 100 markers in each 1 Mb window.")
        nmarkers=length(readdlm(marker_effect_file,',',header=true)[2])
        mapfile = DataFrame(markerID=1:nmarkers,
                            chromosome=fill(1,nmarkers),
                            position=1:10_000:nmarkers*10_000)
        CSV.write("mapfile.temp",mapfile)
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
        for i in 1:length(marker_effects_file)
            #using marker effect files
            output            = readdlm(marker_effects_file[i],',',header=true)[1]
            nsamples,nMarkers = size(output)
            nWindows          = length(window_size_nSNPs)
            winVarProps       = zeros(nsamples,nWindows)
            winVar            = zeros(nsamples,nWindows)
            #window_mrk_start ID and window_mrk_end ID are not provided now
            X = (typeof(mme) <: Array ? mme : mme.M[1].output_genotypes)
            @showprogress "running GWAS..." for i=1:nsamples
                α = output[i,:]
                genVar = var(X*α)
                for winj = 1:length(window_column_start)
                  wStart = window_column_start[winj]
                  wEnd   = window_column_end[winj]
                  BV_winj= X[:,wStart:wEnd]*α[wStart:wEnd]
                  var_winj = var(BV_winj)
                  winVar[i,winj]      = var_winj
                  winVarProps[i,winj] = var_winj/genVar
                end
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
            nsamples,nMarkers = size(output1)
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
             output_winVarProps=false)
     GWAS(mme,map_file,marker_effects_file,
                  header=header,
                  window_size=window_size,
                  threshold=threshold,
                  sliding_window = sliding_window,
                  output_winVarProps=output_winVarProps)
end


export GWAS
