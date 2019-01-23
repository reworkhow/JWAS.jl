function traceplot(file;header=true,separator=',',backend="pyplot",nplots=4)
    #catch errors when no backends are installed
    if backend == "pyplot"
        pyplot()
    elseif backend == "plotly"
        plotly()
    elseif backend == "gr"
        gr()
    end
    mychain,mylabel=readdlm(file,separator,header=header)
    if nplots > size(mychain,2)
        nplots=size(mychain,2)
    end
    mychain = mychain[:,1:nplots]
    mylabel = mylabel[1:nplots]

    steps = 1:size(mychain,1)
    plot(mychain, layout=(nplots,1),title= reshape(mylabel,1,length(mylabel)),
         label="",title_location=:left)
    plot!(cumsum(mychain,dims=1)./steps,layout=(nplots,1),label="",color=:red)
    #title!("trace plot for "*split(file,['/','\\','.'])[end-1])

    # upscale = 8 #8x upscaling in resolution
    # fntsm = Plots.font("sans-serif", 10.0*upscale)
    # fntlg = Plots.font("sans-serif", 14.0*upscale)
    # default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
    # default(size=(800*upscale,600*upscale)) #Plot canvas size
    # default(dpi=300) #Only for PyPlot - presently broken
end

function PSRF(files::AbstractString...;header=true)
    nchain = length(files)

    sample_mean     = []
    sample_variance = []
    chainlength     = 0

    for i in files
        if header == true
            mysample,myheader = readdlm(i,header=header)
        else
            mysample = readdlm(i,header=header)
        end
        push!(sample_mean, mean(mysample))
        push!(sample_variance, std(mysample)^2)
        chainlength = length(mysample)
    end

    M, N = nchain, chainlength
    B = N / (M - 1) * sum((sample_mean .- mean(sample_mean)).^2) # between chain variance
    W = mean(sample_variance)                                    # with-in chain variance

    # unbiased estimator of the marginal posterior variance
    V = (N-1)/N * W + (M + 1)/N*M * B

    return V/W
end
