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
