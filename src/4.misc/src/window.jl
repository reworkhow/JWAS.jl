function GWAS_window(marker_file,mme;header=true,window_size=100)
    if header==true
        output=readdlm(marker_file,header=true)[1]
    else
        output=readdlm(marker_file)
    end

    nsamples,nMarkers=size(output)
    nWindows    = round(Int64,nMarkers/window_size)
    winVarProps = zeros(nsamples,nWindows)
    X           = mme.M.genotypes


    for i=1:nsamples
        α = output[i,:]
        genVar = var(X*α)

        wEnd = 0
        for win=1:nWindows
          wStart = wEnd + 1
          wEnd  += window_size
          wEnd   = (wEnd > nMarkers) ? nMarkers:wEnd
          winVarProps[i,win] = var(X[:,wStart:wEnd]*α[wStart:wEnd])/genVar
        end
    end
    winVarProps
end

export GWAS_window
