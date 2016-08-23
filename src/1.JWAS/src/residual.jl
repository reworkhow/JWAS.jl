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
#make ResVar, dictionary for Rinv
function mkRi(mme::MME,df::DataFrame)
    resVar = ResVar(mme.R,Dict())
    tstMsng = !isna(df[mme.lhsVec[1]])
    #find all missing patterns in data
    for i=2:size(mme.lhsVec,1)
        tstMsng = [tstMsng !isna(df[mme.lhsVec[i]])]
    end
    mme.missingPattern = tstMsng
    n    = size(tstMsng,2)
    nObs = size(tstMsng,1)
    ii = Array(Int64,nObs*n^2)
    jj = Array(Int64,nObs*n^2)
    vv = Array(Float64,nObs*n^2)
    pos = 1
    for i=1:size(tstMsng,1)
        sel = reshape(tstMsng[i,:],n)
        Ri  = getRi(resVar,sel)
        for ti=1:n
            tii = (ti-1)*nObs + i
            for tj=1:n
                tjj = (tj-1)*nObs + i
                ii[pos] = tii
                jj[pos] = tjj
                vv[pos] = Ri[ti,tj]
                pos += 1
            end
        end
    end
    mme.resVar = resVar
    return sparse(ii,jj,vv)
end

function sampleMissingResiduals(mme,resVec)
    msngPtrn = mme.missingPattern
    n,k = size(msngPtrn)
    yIndex = collect(0:k-1)*n
    allTrue = fill(true,k)
    for i=1:n
        notMsng = reshape(msngPtrn[i,:,],k)
        if (notMsng!=allTrue)
            msng    = !notMsng
            nMsng   = sum(msng)
            resi    = resVec[yIndex+i][notMsng]
            Ri      = mme.resVar.RiDict[notMsng][notMsng,notMsng]
            Rc      = mme.R[msng,notMsng]
            L       = chol(mme.R[msng,msng] - Rc*Ri*Rc')'
            resVec[(yIndex+i)[msng]] = Rc*Ri*resi + L*randn(nMsng)
            #resVec[yIndex+i][msng] = Rc*Ri*resi + L*randn(nMsng) this does not work!
        end
    end
end
