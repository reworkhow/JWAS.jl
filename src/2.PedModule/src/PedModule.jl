module PedModule

using DataFrames,CSV, DelimitedFiles
using SparseArrays
using ProgressMeter

export get_pedigree
export get_info

"""
    get_pedigree(pedfile::AbstractString;header=false,separator=',',missingstrings=["0"])
* Get pedigree informtion from a pedigree file with **header** (defaulting to `false`)
  , **separator** (defaulting to `,`) and missing values (defaulting to ["0"])
* Pedigree file format:

```
a,0,0
c,a,b
d,a,c
```
"""
function get_pedigree(pedfile::Union{AbstractString,DataFrames.DataFrame};header=false,separator=',',missingstrings=["0"])
    if typeof(pedfile) <: AbstractString
        printstyled("The delimiter in ",split(pedfile,['/','\\'])[end]," is \'",separator,"\'.\n",bold=false,color=:green)
        df  = CSV.read(pedfile,DataFrame,types=[String,String,String],
                        delim=separator,header=header,missingstrings=missingstrings)
    elseif typeof(pedfile) == DataFrames.DataFrame
        df  = pedfile
    else
        error("Please provide a file path or a dataframe.")
    end
    df[!,1]=strip.(string.(df[!,1]))
    df[!,2]=strip.(string.(df[!,2]))
    df[!,3]=strip.(string.(df[!,3]))
    ped = Pedigree(1,Dict{AbstractString,PedNode}(),
                     Dict{Int64, Float64}(),
                     Set(),Set(),Set(),Set(),Array{String,1}())
    fillMap!(ped,df)
    @showprogress "coding pedigree... " for id in keys(ped.idMap)
     code!(ped,id)
    end
    @showprogress "calculating inbreeding... " for id in keys(ped.idMap)
      calcInbreeding!(ped,id)
    end

    ped.IDs=getIDs(ped)

    get_info(ped)
    writedlm("IDs_for_individuals_with_pedigree.txt",ped.IDs)

    return ped
end

mutable struct PedNode
    seqID::Int64
    sire::String
    dam::String
    f::Float64
end

mutable struct Pedigree
    currentID::Int64
    idMap::Dict{AbstractString,PedNode}
    aij::Dict{Int64, Float64}
    setNG::Set
    setG::Set
    setG_core::Set
    setG_notcore::Set
    #individuals IDs (in order of numerator relationship matrix A)
    IDs::Array{String,1}
end

function code!(ped::Pedigree,id::AbstractString)
# The idea for this function came from a perl script by Bernt Guldbrandtsen
    if ped.idMap[id].seqID!=0
        return
    end
    sireID = ped.idMap[id].sire
    damID  = ped.idMap[id].dam
    if sireID!="missing" && ped.idMap[sireID].seqID==0
        code!(ped,sireID)
    end
    if damID!="missing" && ped.idMap[damID].seqID==0
        code!(ped,damID)
    end
    ped.idMap[id].seqID = ped.currentID
    ped.currentID += 1
end

function fillMap!(ped::Pedigree,df)
    #Warning: indexing with colon as row will create a copy in the future
    #use df[col_inds] to get the columns without copying
    n = size(df,1)
    for i in df[!,2] #same to df[:,2] in deprecated CSV
        if i!="missing" && !haskey(ped.idMap,i)          # skip 0 and if already done
            ped.idMap[i]=PedNode(0,"missing","missing",-1.0)
        end
    end
    for i in df[!,3]
        if i!="missing" && !haskey(ped.idMap,i)         # make an entry for all dams
            ped.idMap[i]=PedNode(0,"missing","missing",-1.0)
        end
    end
    j=1
    for i in df[!,1]
        ped.idMap[i]=PedNode(0,df[j,2],df[j,3],-1.0)
        j+=1
    end
end

function calcAddRel!(ped::Pedigree,id1::AbstractString,id2::AbstractString)
    #@printf "calcRel between %s and %s \n" id1 id2
    if id1=="missing" || id2=="missing"           # zero
        return 0.0
    end
    old,yng = ped.idMap[id1].seqID < ped.idMap[id2].seqID ? (id1,id2) : (id2,id1)
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

    aOldDamYoung  = (old=="missing" || damOfYng =="missing") ? 0.0 : calcAddRel!(ped,old,damOfYng)
    aOldSireYoung = (old=="missing" || sireOfYng=="missing") ? 0.0 : calcAddRel!(ped,old,sireOfYng)
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
    if (sireID=="missing" || damID=="missing" )
        ped.idMap[id].f = 0.0
    else
        ped.idMap[id].f = 0.5*calcAddRel!(ped,sireID,damID)
    end
end

function AInverse(ped::Pedigree)
    ii,jj,vv = HAi(ped)
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
        sirePos = sire=="missing" ? 0 : ped.idMap[sire].seqID
        damPos  = dam =="missing" ? 0 : ped.idMap[dam ].seqID
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

function AInverseSlow(ped::Pedigree)
    n = ped.currentID - 1
    Ai = spzeros(n,n)
    pos  = Int64[0,0,0]
    q    = [0.5,0.5,1.0]
    for ind in keys(ped.idMap)
        sire = ped.idMap[ind].sire
        dam  = ped.idMap[ind].dam
        pos[1] = sire=="missing" ? 0 : ped.idMap[sire].seqID
        pos[2] = dam =="missing" ? 0 : ped.idMap[dam ].seqID
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

function getIDs(ped::Pedigree)
    n = length(ped.idMap)
    ids = Array{String}(undef,n)
    for i in ped.idMap
      ids[i[2].seqID] = i[1]
    end
    return ids
end

function getInbreeding(ped::Pedigree)
    n = length(ped.idMap)
    inbreeding = Array{AbstractFloat}(undef,n)
    for i in ped.idMap
      inbreeding[i[2].seqID] = i[2].f
    end
    return inbreeding
end

"""
    get_info(pedigree::Pedigree;Ai=false)
* Print summary informtion from a pedigree object including number of individulas, sires.
  dams and founders. Return individual IDs, inverse of numerator relationship matrix,
  and inbreeding coefficients if **Ai**=`true`.

"""
function get_info(pedigree::Pedigree;Ai=false)
    println("Pedigree information:")
    println("#individuals: ",length(pedigree.idMap))
    sires  = [pednode.sire for pednode in values(pedigree.idMap)]
    dams = [pednode.dam for pednode in values(pedigree.idMap)]
    println("#sires:       ",sum(unique(sires).!="missing"))
    println("#dams:        ",sum(unique(dams).!="missing"))
    println("#founders:    ",length(pedigree.idMap)-sum((sires .!= "missing") .* (dams .!= "missing")))

    if Ai == true
        Ai   = PedModule.AInverse(pedigree)
        IDs = PedModule.getIDs(pedigree)
        inbreeding = PedModule.getInbreeding(pedigree)
        println("Get individual IDs, inverse of numerator relationship matrix, and inbreeding coefficients.")
        return IDs,Ai,inbreeding
    end
end

include("forSSBR.jl")

end # of PedModule
