#GBLUP functions
function GBLUP_setup(Mi::Genotypes)
    Mi.genotypes  = Mi.genotypes ./ sqrt.(2*Mi.alleleFreq.*(1 .- Mi.alleleFreq))
    G       = (Mi.genotypes*Mi.genotypes'+ I*0.00001)/Mi.nMarkers
    eigenG  = eigen(G)
    L       = eigenG.vectors
    D       = eigenG.values
    # α is pseudo marker effects of length nobs (starting values = L'(starting value for BV)
    Mi.nMarkers= Mi.nObs
    #reset parameters in output
    M2   = Mi.output_genotypes ./ sqrt.(2*Mi.alleleFreq.*(1 .- Mi.alleleFreq))
    M2Mt = M2*Mi.genotypes'/Mi.nMarkers
    Mi.output_genotypes = M2Mt*L*Diagonal(1 ./D)
    #reset parameter in mme.M
    Mi.G         = Mi.genetic_variance
    Mi.scale     = Mi.G*(Mi.df-2)/Mi.df
    Mi.markerID  = string.(1:Mi.nObs) #pseudo markers of length=nObs
    Mi.genotypes = L
    Mi.α         = L'Mi.α
    return D
end

function GBLUP!(Mi::Genotypes,ycorr,vare,Rinv,D)
    ycorr       = ycorr + Mi.genotypes*Mi.α
    lhs         = Rinv .+ vare./(Mi.G*D)
    mean1       = Mi.genotypes'*(Rinv.*ycorr)./lhs
    Mi.α        = mean1 + randn(Mi.nObs).*sqrt.(vare./lhs)
    ycorr[:]    = ycorr - Mi.genotypes*Mi.α
end
