"""
    GWAS(marker_effects_file,model;header=true,window_size=100,threshold=0.001)

run genomic window-based GWAS without a map file

* MCMC samples of marker effects are stored in **marker_effects_file** with delimiter ','.
* **window_size** is either a constant (identical number of markers in each window) or an array of number of markers in each window
* **model** is either the model::MME used in analysis or the genotypic covariate matrix M::Array
* File format:

"""
function GWAS(marker_effect_file,mme;header=true,window_size=100,threshold=0.001,output_winVarProps=false)
    println("This function is deprecated. Please check ?GWAS")

    if typeof(mme) <: Array
        nMarkers = size(mme,2)
    else
        nMarkers = size(mme.M.markerID,1)
    end
    if length(window_size)==1
        nWindows    = ceil(Int64,nMarkers/window_size)
        window_size = fill(window_size,nWindows)
    end

    winVarProps, window_mrk_start, window_mrk_end = getWinVarProps(marker_effect_file,mme;header=header,window_size=window_size,threshold=threshold)

    winVarProps[isnan.(winVarProps)] .= 0.0 #replace NaN caused by situations no markers are included in the model
    WPPA, prop_genvar = vec(mean(winVarProps .> threshold,dims=1)), vec(mean(winVarProps,dims=1))
    prop_genvar = round.(prop_genvar*100,digits=2)

    srtIndx = sortperm(WPPA,rev=true)
    out = DataFrame(wStart = window_mrk_start[srtIndx],
                    wEnd   = window_mrk_end[srtIndx],
                    wSize  = window_size[srtIndx],
                    prGenVar = prop_genvar[srtIndx],
                    WPPA     = WPPA[srtIndx],
                    PPA_t  = cumsum(WPPA[srtIndx]) ./ (1:length(WPPA))
                   )
    if output_winVarProps == false
        return out
    else
        return out,winVarProps
    end
end

function getWinVarProps(marker_effect_file,mme;header=true,window_size=false,threshold=0.001)
    if header==true
        output=readdlm(marker_effect_file,',',header=true)[1]
    else
        output=readdlm(marker_effect_file,',')
    end

    nsamples,nMarkers=size(output)
    nWindows    = length(window_size)

    winVarProps = zeros(nsamples,nWindows)
    if typeof(mme) <: Array
        X = mme
        window_mrk_start = Array{Int64,1}()
        window_mrk_end   = Array{Int64,1}()
    else
        X = mme.output_genotypes
        window_mrk_start = Array{String,1}()
        window_mrk_end   = Array{String,1}()
    end

    wEnd    = 0
    @showprogress  for i in window_size
        wStart = wEnd + 1
        wEnd  += i
        wEnd   = (wEnd > nMarkers) ? nMarkers : wEnd
        if typeof(mme) <: Array
            push!(window_mrk_start,wStart)
            push!(window_mrk_end  ,wEnd)
        else
            push!(window_mrk_start,mme.M.markerID[wStart])
            push!(window_mrk_end  ,mme.M.markerID[wEnd])
        end
    end

    for i=1:nsamples
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
    return winVarProps, window_mrk_start, window_mrk_end
end
