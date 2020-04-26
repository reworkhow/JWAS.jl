#GBLUP functions
function GBLUP_setup(Mi::Genotypes) #for both single-trait and multi-trait analysis
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
    if Mi.ntraits == 1
        Mi.scale = Mi.G*(Mi.df-Mi.ntraits-1)/Mi.df
    else
        Mi.scale = Mi.G*(Mi.df-Mi.ntraits-1)
    end
    Mi.markerID  = string.(1:Mi.nObs) #pseudo markers of length=nObs
    Mi.genotypes = L
    Mi.α         = L'Mi.α
    Mi.D         = D
end

function GBLUP!(Mi::Genotypes,ycorr,vare,Rinv)
    ycorr       = ycorr + Mi.genotypes*Mi.α
    lhs         = Rinv .+ vare./(Mi.G*Mi.D)
    mean1       = Mi.genotypes'*(Rinv.*ycorr)./lhs
    Mi.α        = mean1 + randn(Mi.nObs).*sqrt.(vare./lhs)
    ycorr[:]    = ycorr - Mi.genotypes*Mi.α
end


function MTGBLUP!(Mi::Genotypes,ycorr_array,ycorr,vare,Rinv)
    iR0      = inv(vare)
    iGM      = inv(Mi.G)
    for trait = 1:Mi.ntraits
        ycorr_array[trait][:] = ycorr_array[trait] + Mi.genotypes*Mi.α[trait]
    end
    lhs    = [iR0*Rinv[i] + iGM/Mi.D[i] for i=1:length(Mi.D)]
    RHS    = (Mi.genotypes'Diagonal(Rinv))*reshape(ycorr,Mi.nObs,Mi.ntraits)*iR0 #size nmarkers (=nObs) * nTraits
    rhs    = [RHS[i,:] for i in 1:size(RHS,1)]  #not column major
    Σα     = Symmetric.(inv.(lhs))              #nTrait X nTrait
    μα     = Σα.*rhs
    αs     = rand.(MvNormal.(μα,Σα))
    for markeri = 1:Mi.nMarkers
        for traiti = 1:Mi.ntraits
            Mi.α[traiti][markeri] = αs[markeri][traiti]
        end
    end
    for trait = 1:Mi.ntraits
        ycorr_array[trait][:] = ycorr_array[trait] - Mi.genotypes*Mi.α[trait]
    end
end
