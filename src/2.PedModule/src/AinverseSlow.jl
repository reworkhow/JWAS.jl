function AInverseSlow(ped::Pedigree)
    n = ped.currentID - 1
    Ai = spzeros(n,n)
    pos  = Int64[0,0,0]
    q    = [0.5,0.5,1.0]
    for ind in keys(ped.idMap)
        sire = ped.idMap[ind].sire
        dam  = ped.idMap[ind].dam
        pos[1] = sire=="0" ? 0: ped.idMap[sire].seqID
        pos[2] = dam =="0" ? 0: ped.idMap[dam ].seqID
        pos[3] = ped.idMap[ind].seqID
        if pos[1]>0 && pos[2]>0
            q[1] = -0.5
            q[2] = -0.5
            d = 4.0/(2 - ped.idMap[sire].f - ped.idMap[dam].f)
        elseif pos[1]>0
            q[1] = -0.5
            q[2] = 0.0
            d = 4.0/(3 - ped.idMap[sire].f)
        elseif pos[2]>0
            q[1] = 0.0
            q[2] = -0.5
            d = 4.0/(3 - ped.idMap[dam].f)
        else
            q[1] = 0.0
            q[2] = 0.0
            d = 1.0
        end
        for i=1:3
            ii = pos[i]
            if ii>0
                for j=1:3
                    jj = pos[j]
                    if jj>0
                        Ai[ii,jj] += q[i]*q[j]*d
                    end
                end
            end
        end
    end
    return (Ai)
end