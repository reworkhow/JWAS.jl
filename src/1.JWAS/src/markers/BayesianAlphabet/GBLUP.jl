################################################################################
#GBLUP: pseudo markers are used.
#y = μ + a + e with mean(a)=0,var(a)=Gσ²=MM'σ² and G = LDL' <==>
#y = μ + Lα +e with mean(α)=0,var(α)=D*σ² : L orthogonal
#y2hat = cov(a2,a)*inv(var(a))*L*αhat =>
#      = (M2*M')*inv(MM')*L*αhat = (M2*M')*(L(1./D)L')(Lαhat)=(M2*M'*L(1./D))αhat
################################################################################
#GBLUP functions
function GBLUP_setup(Mi::Genotypes) #for both single-trait and multi-trait analysis
    G       = Mi.genotypes
    eigenG  = eigen(G)
    L       = eigenG.vectors
    D       = eigenG.values
    # α is pseudo marker effects of length nobs (starting values = L'(starting value for BV)
    Mi.nMarkers= Mi.nObs
    #reset parameters in output
    if Mi.isGRM #if Genomic relationship matrix is provided,
        M2Mt  = Mi.output_genotypes
    else                      #calculate the relationship matrix from the genotype covariate matrix
        M2   = Mi.output_genotypes ./ sqrt.(2*Mi.alleleFreq.*(1 .- Mi.alleleFreq))
        M2Mt = M2*Mi.genotypes'/Mi.nMarkers
    end
    Mi.output_genotypes = M2Mt*L*Diagonal(1 ./D)
    #reset parameter in mme.M
    Mi.markerID  = string.(1:Mi.nObs) #pseudo markers of length=nObs
    Mi.genotypes = L
    for traiti = 1:Mi.ntraits
        Mi.α[traiti] = L'Mi.α[traiti]
    end
    Mi.D         = abs.(D) #avoid very small negative values
end

function megaGBLUP!(Mi::Genotypes,wArray,vare,Rinv)
    Threads.@threads for i in 1:length(wArray) #ntraits
        GBLUP!(Mi.genotypes,Mi.α[i],Mi.D,wArray[i],vare[i,i],Mi.G.val[i,i],Rinv,Mi.nObs)
    end
end

function GBLUP!(Mi::Genotypes,ycorr,vare,Rinv) #single-trait
    GBLUP!(Mi.genotypes,Mi.α[1],Mi.D,ycorr,vare,Mi.G.val[1,1],Rinv,Mi.nObs)
end

function GBLUP!(genotypes,α,D,ycorr,vare,vara,Rinv,nObs)
    ycorr[:]    = ycorr + genotypes*α #ycor[:] is needed (ycor causes problems)
    lhs         = Rinv .+ vare./(vara*D)
    mean1       = genotypes'*(Rinv.*ycorr)./lhs
    α[:]        = mean1 + randn(nObs).*sqrt.(vare./lhs)
    ycorr[:]    = ycorr - genotypes*α
end


function MTGBLUP!(Mi::Genotypes,ycorr_array,ycorr,vare,Rinv)
    iR0      = inv(vare)
    iGM      = inv(Mi.G.val)
    for trait = 1:Mi.ntraits
        ycorr_array[trait][:] = ycorr_array[trait] + Mi.genotypes*Mi.α[trait]
    end
    lhs    = [iR0*Rinv[i] + iGM/Mi.D[i] for i=1:length(Mi.D)]
    RHS    = (Mi.genotypes'Diagonal(Rinv))*reshape(ycorr,Mi.nObs,Mi.ntraits)*iR0 #size nmarkers (=nObs) * ntraits
    rhs    = [RHS[i,:] for i in 1:size(RHS,1)]  #not column major
    Σα     = Symmetric.(inv.(lhs))              #ntraits X ntraits
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
