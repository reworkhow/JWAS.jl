#fill up missing phenotype patterns with corresponding inverted residual variance
function getRi(resVar::ResVar,sel::BitArray{1})
    if haskey(resVar.RiDict,sel)
        return resVar.RiDict[sel]
    end
    n = size(resVar.R0,1)
    RZ = zeros(n,n)
    RZ[sel,sel] = inv(resVar.R0[sel,sel])
    resVar.RiDict[sel] = copy(RZ)
    return RZ
end

#make tricky Ri (big) allowing NA in phenotypes and fixed effects
#make ResVar, dictionary for Rinv, no sample Missing residuals
function mkRi(mme::MME,df::DataFrame,Rinv)
    resVar   = ResVar(mme.R.val,Dict())
    tstMsng  = .!ismissing.(Matrix(df[!,mme.lhsVec]))
    if mme.missingPattern == false
        mme.missingPattern = tstMsng
    end
    ntrait = size(tstMsng,2)
    nObs   = size(tstMsng,1)
    ii = Array{Int64}(undef,nObs*ntrait^2)
    jj = Array{Int64}(undef,nObs*ntrait^2)
    vv = Array{AbstractFloat}(undef,nObs*ntrait^2)
    pos = 1
    for i=1:nObs
        sel = tstMsng[i,:]
        Ri  = getRi(resVar,sel)*Rinv[i]
        for ti=1:ntrait
            tii = (ti-1)*nObs + i
            for tj=1:ntrait
                tjj = (tj-1)*nObs + i
                ii[pos] = tii
                jj[pos] = tjj
                vv[pos] = Ri[ti,tj]
                pos += 1
            end
        end
    end
    mme.resVar = resVar
    # This narrows the element type of the vv vector
    vv = map(identity, vv)
    return sparse(ii,jj,vv)
end


#For other methods especially those with marker effects
#Missing residuals (phenotypes) are imputed to enable estimation of residual
#variances and adjusting phenotyp approach to sample marker effects.
function sampleMissingResiduals(mme,resVec)
    msngPtrn = mme.missingPattern
    nobs,ntraits = size(msngPtrn)
    yIndex = collect(0:(ntraits-1))*nobs
    allTrue = fill(true,ntraits)
    mydata = reshape(resVec,nobs,ntraits)
    for (notmissing, _ ) in mme.resVar.RiDict
        if notmissing!=allTrue
            #find all inds with same missing pattern
            mybool = msngPtrn[:,1].==notmissing[1]
            for i in 2:size(msngPtrn,2)
                mybool = mybool .& (msngPtrn[:,i].==notmissing[i])
            end
            #create cov and var matrices for BLP imputation
            Ri  = inv(Symmetric(mme.R.val[notmissing,notmissing]))
            Rc  = mme.R.val[.!notmissing,notmissing]
            L   = (cholesky(Symmetric(mme.R.val[.!notmissing,.!notmissing] - Rc*Ri*Rc')).U)#scalar cholesky
            #imputation
            mydata[mybool,.!notmissing]= mydata[mybool,notmissing]*Ri*Rc' + randn(sum(mybool),sum(.!notmissing))*L
        end
    end
    return reshape(mydata,nobs*ntraits)
end

#SLOW DUE TO ACCESSING A SUBSET OF A MATRIX MANY TIMES
# function sampleMissingResiduals(mme,resVec)
#     msngPtrn = mme.missingPattern
#     nobs,ntraits = size(msngPtrn)
#     yIndex = collect(0:(ntraits-1))*nobs
#     allTrue = fill(true,ntraits)
#     for i=1:nobs
#         notMsng = msngPtrn[i,:]
#         if (notMsng!=allTrue)
#             msng    = .!notMsng
#             resi    = resVec[yIndex .+ i][notMsng]
#             Ri      = mme.resVar.RiDict[notMsng][notMsng,notMsng]
#             Rc      = mme.R.val[msng,notMsng]
#             L       = (cholesky(Symmetric(mme.R.val[msng,msng] - Rc*Ri*Rc')).U)'
#             resVec[(yIndex .+ i)[msng]] = Rc*Ri*resi + L*randn(sum(msng))
#             #resVec[yIndex+i][msng] = Rc*Ri*resi + L*randn(nMsng) this does not work!
#         end
#     end
# end
