module PedModule

using DataFrames

type PedNode
    seqID::Int64
    sire::ASCIIString
    dam::ASCIIString
    f::Float64
end

type Pedigree
    currentID::Int64
    idMap::Dict{AbstractString,PedNode}
    aij::Dict{Int64, Float64}
    setNG::Set
    setG::Set
    setG_core::Set
    setG_notcore::Set
end

function code!(ped::Pedigree,id::AbstractString)
# The idea for this function came from a perl script by Bernt Guldbrandtsen
    if ped.idMap[id].seqID!=0
        return
    end
    sireID = ped.idMap[id].sire
    damID  = ped.idMap[id].dam
    if sireID!="0" && ped.idMap[sireID].seqID==0
        code!(ped,sireID)
    end
    if damID!="0" && ped.idMap[damID].seqID==0
        code!(ped,damID)
    end
    ped.idMap[id].seqID = ped.currentID
    ped.currentID += 1
end

function fillMap!(ped::Pedigree,df)
    n = size(df,1)
    for i in df[:,2]
        if i!="0" && !haskey(ped.idMap,i)          # skip 0 and if already done
            ped.idMap[i]=PedNode(0,"0","0",-1.0)
        end
    end
    for i in df[:,3]
        if i!="0" && !haskey(ped.idMap,i)         # make an entry for all dams
            ped.idMap[i]=PedNode(0,"0","0",-1.0)
        end
    end
    j=1
    for i in df[:,1]
        ped.idMap[i]=PedNode(0,df[j,2],df[j,3],-1.0)
        j+=1
    end
end

function calcAddRel!(ped::Pedigree,id1::AbstractString,id2::AbstractString)
    #@printf "calcRel between %s and %s \n" id1 id2
    if id1=="0" || id2=="0"           # zero
        return 0.0
    end
    old,yng = ped.idMap[id1].seqID<ped.idMap[id2].seqID ? (id1,id2):(id2,id1)
    oldID = ped.idMap[old].seqID
    yngID = ped.idMap[yng].seqID
    
    n = yngID - 1
    aijKey = n*(n+1)/2 + oldID 
    if haskey(ped.aij,aijKey)     
        return ped.aij[aijKey]
    end       
    
    sireOfYng = ped.idMap[yng].sire
    damOfYng  = ped.idMap[yng].dam
    
    if old==yng                       # aii
        #aii = 1.0 + calcInbreeding!(ped,old)
        aii = 1.0 + 0.5*calcAddRel!(ped,sireOfYng,damOfYng)
        ped.aij[aijKey] = aii
        return (aii)
    end
    
    aOldDamYoung  = (old=="0" || damOfYng =="0")? 0.0:calcAddRel!(ped,old,damOfYng)
    aOldSireYoung = (old=="0" || sireOfYng=="0")? 0.0:calcAddRel!(ped,old,sireOfYng)
    aijVal = 0.5*(aOldSireYoung + aOldDamYoung)
    ped.aij[aijKey] = aijVal

    #aij = 0.5*(calcAddRel!(ped,old,sireOfYng) + calcAddRel!(ped,old,damOfYng))
    #ped.aij[yngID,oldID] = aij
    #ped.aij[oldID,yngID] = 1.0
    return aijVal
end

function calcInbreeding!(ped::Pedigree,id::AbstractString)
    #@printf "calcInbreeding for: %s \n" id
    if ped.idMap[id].f > -1.0
        return ped.idMap[id].f
    end
    sireID = ped.idMap[id].sire
    damID  = ped.idMap[id].dam
    if (sireID=="0" || damID=="0" )
        ped.idMap[id].f = 0.0
    else
        ped.idMap[id].f = 0.5*calcAddRel!(ped,sireID,damID)
    end
end

function AInverse(ped::Pedigree)
    ii,jj,vv = HAi(ped)#cholesky??
    hAi      = sparse(ii,jj,vv)
    Ai       = hAi'hAi
    return Ai
end

function HAi(ped::Pedigree)
    ii = Int64[]
    jj = Int64[]
    vv = Float64[]
    for ind in keys(ped.idMap)
        sire = ped.idMap[ind].sire
        dam  = ped.idMap[ind].dam
        sirePos = sire=="0" ? 0: ped.idMap[sire].seqID
        damPos  = dam =="0" ? 0: ped.idMap[dam ].seqID
        myPos   = ped.idMap[ind].seqID
        if sirePos>0 && damPos>0
            d = sqrt(4.0/(2 - ped.idMap[sire].f - ped.idMap[dam].f))
            push!(ii,myPos)
            push!(jj,sirePos)
            push!(vv,-0.5*d)
            push!(ii,myPos)
            push!(jj,damPos)
            push!(vv,-0.5*d)
            push!(ii,myPos)
            push!(jj,myPos)
            push!(vv,d)
         elseif sirePos>0
            d = sqrt(4.0/(3 - ped.idMap[sire].f))
            push!(ii,myPos)
            push!(jj,sirePos)
            push!(vv,-0.5*d)
            push!(ii,myPos)
            push!(jj,myPos)
            push!(vv,d)
          elseif damPos>0
            d = sqrt(4.0/(3 - ped.idMap[dam].f))
            push!(ii,myPos)
            push!(jj,damPos)
            push!(vv,-0.5*d)
            push!(ii,myPos)
            push!(jj,myPos)
            push!(vv,d)
        else
            d = 1.0
            push!(ii,myPos)
            push!(jj,myPos)
            push!(vv,d)
        end
    end
    return (ii,jj,vv)
end

function  mkPed(pedFile::AbstractString;header=false,separator=' ')
    #dataframes string conflits with AbstractString in julia(fixed)
    #df = readtable(pedFile,separator = ' ',header=false)

    df  = readtable(pedFile,eltypes=[UTF8String,UTF8String,UTF8String],
                            separator=separator,header=header)
    ped = Pedigree(1,Dict{AbstractString,PedNode}(),
                     Dict{Int64, Float64}(),
                     Set(),Set(),Set(),Set())

    fillMap!(ped,df)
    
    for id in keys(ped.idMap)
     code!(ped,id)
    end
    
    for id in keys(ped.idMap)
      calcInbreeding!(ped,id)
    end
    
    return ped
end

function getIDs(ped::Pedigree)
    n = length(ped.idMap)
    ids = Array(ASCIIString,n)
    for i in ped.idMap
      ids[i[2].seqID] = i[1]
    end
    return ids
end

include("forSSBR.jl")

end # of PedModule
