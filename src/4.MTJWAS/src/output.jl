function GWAS_window(file,mme::MME;windowSize=100)
    output=readdlm(file)
    nsamples,nMarkers=size(output)
    nWindows    = round(Int64,nMarkers/windowSize)
    winVarProps = zeros(nsamples,nWindows)
    X           = mme.M.X
    
    
    for i=1:nsamples
        α = output[i,:]'
        genVar = var(X*α)
            
        wEnd = 0
        for win=1:nWindows
          wStart = wEnd + 1
          wEnd  += windowSize
          wEnd   = (wEnd > nMarkers) ? nMarkers:wEnd
          winVarProps[i,win] = var(X[:,wStart:wEnd]*α[wStart:wEnd])#/genVar
        end
    end
    winVarProps
end