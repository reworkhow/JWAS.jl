"""
    GWAS(marker_effects_file;header=false)

Compute the model frequency for each marker (the probability the marker is included in the model) using samples of marker effects stored in **marker_effects_file**.
"""
function GWAS(marker_effects_file;header=true)
    file = marker_effects_file
    if header==true
        samples,markerID =readdlm(file,',',header=true)
    else
        samples=readdlm(file,',')
    end

    modelfrequency = vec(mean(samples .!= 0.0,dims=1))
    if header==true
        out = [vec(markerID) modelfrequency]
    else
        out = modelfrequency
    end
    return out
end

"""
    GWAS(marker_effects_file,model;header=true,window_size="1 Mb",threshold=0.001)

run genomic window-based GWAS without marker locations

* MCMC samples of marker effects are stored in **marker_effects_file** with delimiter ','.
* **window_size** is either a constant (identical number of markers in each window) or an array of number of markers in each window
* **model** is either the model::MME used in analysis or the genotypic cavariate matrix M::Array
* File format:

"""
function GWAS(marker_effect_file,mme;header=true,window_size=100,threshold=0.001)
    println("Compute the posterior probability of association of the genomic window that explains more than ",threshold," of the total genetic variance")

    if header==true
        output=readdlm(marker_effect_file,',',header=true)[1]
    else
        output=readdlm(marker_effect_file,',')
    end

    nsamples,nMarkers=size(output)
    if length(window_size)==1
        nWindows    = ceil(Int64,nMarkers/window_size)
        window_size = fill(window_size,nWindows)
    else
        nWindows    = length(window_size)
    end

    winVarProps = zeros(nsamples,nWindows)
    if typeof(mme) <: Array
        X = mme
    else
        X = mme.output_genotypes
    end

    @showprogress for i=1:nsamples
        α = output[i,:]
        genVar = var(X*α)

        wEnd    = 0
        windowi = 1
        for win=1:nWindows
          wStart = wEnd + 1
          wEnd  += window_size[windowi]
          wEnd   = (wEnd > nMarkers) ? nMarkers : wEnd
          winVarProps[i,win] = var(X[:,wStart:wEnd]*α[wStart:wEnd])/genVar

          windowi +=1
        end
    end
    return winVarProps
end


#Chen, C., Steibel, J. P., & Tempelman, R. J. (2017). Genome-Wide Association
#Analyses Based on Broadly Different Specifications for Prior Distributions,
#Genomic Windows, and Estimation Methods. Genetics, 206(4), 1791–1806.
#Also, as per Moser et al. (2015), two different within-chromosome starting positions
#(starting at location 0 or 0.25 Mb for window sizes 0.5, starting at 0 or 0.5 Mb
#location for window sizes 1 Mb, starting at 0 or 1 Mb location for window sizes 2 Mb,
#and starting at 0 or 1.5 Mb location for window sizes 3 Mb) for each chromosome were
#chosen to partly counteract the chance effect of different LD patterns being associated
#with nonoverlapping windows. Fi- nally, adaptive window sizes based on clustering SNP
#by LD r2 were also determined using the BALD R package (Dehman and Neuvial 2015),
#using the procedure described by Dehman et al. (2015).
"""
    GWAS(marker_effects_file,map_file,model;header=true,window_size="1 Mb",threshold=0.001)

run genomic window-based GWAS

* MCMC samples of marker effects are stored in **marker_effects_file** with delimiter ','.
* **model** is either the model::MME used in analysis or the genotypic cavariate matrix M::Array
* **map_file** has the (sorted) marker position information with delimiter ','.
* File format:

```
markerID,chromosome,position
m1,1,16977
m2,1,434311
m3,1,1025513
m4,2,70350
m5,2,101135
```


"""
function GWAS(marker_effects_file,map_file,mme;header=true,window_size="1 Mb",threshold=0.001,output_winVarProps=false)

    if window_size=="1 Mb"
        window_size_Mb = 1_000_000
    else
        window_size_Mb = map(Int64,parse(Float64,split(window_size)[1])*1_000_000)
    end

    if header == true
        mapfile = readdlm(map_file,',',header=true)[1]
    else
        mapfile = readdlm(map_file,',')
    end
    chr     = map(string,mapfile[:,2])
    pos     = map(Int64,mapfile[:,3])

    window_size = Array{Int64,1}() #save number of markers in ith window for all windows
    window_chr  = Array{String,1}()
    window_pos_start = Array{Int64,1}()
    window_pos_end   = Array{Int64,1}()
    window_snp_start = Array{Int64,1}()
    window_snp_end   = Array{Int64,1}()

    for i in unique(chr)
      pos_on_chri     = pos[chr.== i] #assume chr and pos are sorted
      nwindow_on_chri = ceil(Int64,pos_on_chri[end]/window_size_Mb)

      for j in 1: nwindow_on_chri
        thisstart= window_size_Mb*(j-1)
        thisend  = window_size_Mb*j
        push!(window_chr,i)
        push!(window_pos_start,thisstart)
        push!(window_pos_end,thisend)
        snps_window = thisstart .< pos_on_chri .<= thisend
        push!(window_size,sum(snps_window))
        if sum(snps_window)!=0
            push!(window_snp_start,pos_on_chri[findfirst(snps_window)])
            push!(window_snp_end,pos_on_chri[findlast(snps_window)])
        else
            push!(window_snp_start,0)
            push!(window_snp_end,0)
        end
      end
    end
    winVarProps = GWAS(marker_effects_file,mme,header=header,window_size=window_size,threshold=threshold)
    winVarProps[isnan.(winVarProps)] .= 0.0 #replace NaN caused by situations no markers are included in the model
    WPPA, prop_genvar = vec(mean(winVarProps .> threshold,dims=1)), vec(mean(winVarProps,dims=1))
    prop_genvar = round.(prop_genvar*100,digits=2)
    #bug in Julia, vcat too long
    #out  = [["window";1:length(WPPA)] ["chr"; window_chr] ["start"; window_pos_start] ["end"; window_pos_end]["WPPA"; WPPA]]
    out1 = [["window";1:length(WPPA)] ["chr"; window_chr]]
    out2 = [["start"; window_pos_start] ["end"; window_pos_end]]
    out3 = [["start_SNP"; window_snp_start] ["end_SNP"; window_snp_end] ["#SNPS"; window_size]]
    out4 = [["%genetic variance";prop_genvar] ["WPPA"; WPPA]]
    out  = [out1 out2 out3 out4]
    if output_winVarProps == false
        return out
    else
        return out,winVarProps
    end
end

export GWAS
