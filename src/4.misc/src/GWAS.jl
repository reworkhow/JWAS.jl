#return model frequencies
function GWAS(marker_effects_file;header=true)
    file = marker_effects_file
    if header==true
        samples,markerID =readdlm(file,header=true)[1]
    else
        samples=readdlm(file)
    end

    modelfrequency = vec(mean(samples .!= 0.0,1))
    if header==true
        out = [markerID modelfrequency]
    else
        out = modelfrequency
    end
    return out
end

#window size
#same number of SNPs in each window
#an array of number of SNPs in each window

function GWAS(marker_file,mme;header=true,window_size=100,threshold=0.001)
    println("Compute posterior probability that window explains more than ",threshold," of genetic variance")

    if header==true
        output=readdlm(marker_file,header=true)[1]
    else
        output=readdlm(marker_file)
    end

    nsamples,nMarkers=size(output)
    if length(window_size)==1
        nWindows    = ceil(Int64,nMarkers/window_size)
        window_size = fill(window_size,nWindows)
    else
        nWindows    = length(window_size)
    end

    winVarProps = zeros(nsamples,nWindows)
    X           = mme.M.genotypes

    @showprogress for i=1:nsamples
        α = output[i,:]
        genVar = var(X*α)

        wEnd    = 0
        windowi = 1
        for win=1:nWindows
          wStart = wEnd + 1
          wEnd  += window_size[windowi]
          wEnd   = (wEnd > nMarkers) ? nMarkers:wEnd
          winVarProps[i,win] = var(X[:,wStart:wEnd]*α[wStart:wEnd])/genVar

          windowi +=1
        end
    end
    #return(vec(mean(winVarProps .> threshold,1)), mean(winVarProps,1))
    return vec(mean(winVarProps .> threshold,1))
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

function GWAS(marker_effects_file,map_file,mme;header=true,window_size="1 Mb",threshold=0.001)

    if window_size=="1 Mb"
        window_size_Mb = 1_000_000
    else
        window_size_Mb = map(Int64,parse(Float64,split(window_size)[1])*1_000_000)
    end

    mapfile = readdlm(map_file)
    chr     = map(Int64,mapfile[:,2])
    pos     = map(Int64,mapfile[:,3])

    window_size = Array{Int64,1}() #save number of markers in ith window for all windows
    for i in 1:maximum(chr)
      pos_on_chri     = pos[chr.==i]
      nwindow_on_chri = ceil(Int64,pos_on_chri[end]/window_size_Mb)

      for j in 1: nwindow_on_chri
        push!(window_size,sum(window_size_Mb*(j-1) .< pos_on_chri .< window_size_Mb*j))
      end
    end
    return GWAS(marker_effects_file,mme,header=header,window_size=window_size,threshold=threshold)
end



export GWAS
