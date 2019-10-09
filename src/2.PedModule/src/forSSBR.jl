function genoSet!(genoID_file::AbstractString,ped::Pedigree)
    df = readtable(genoID_file, eltypes=[String], separator = ' ',header=false)
	df = vec(map(string,convert(Array,df)))
    for i in df
        push!(ped.setG,i)
	end
    all = Set()
	for i in keys(ped.idMap)
        push!(all,i)
	end
    ped.setNG = setdiff(all,ped.setG)
	j = 1
	for i in ped.setNG
		ped.idMap[i].seqID = j
		j += 1
	end
	numberNonGeno = j - 1
	for i in ped.setG
		ped.idMap[i].seqID = j
		j += 1
	end

	ped.IDs=getIDs(ped) #order of IDs is changed

	return (numberNonGeno)
end

function genoSet!(genoID::Array{AbstractString,1},ped::Pedigree)
    if !issubset(genoID,map(string,collect(keys(ped.idMap))))
        error("Not all genotyped individuals are in the pedigree! Please only use genotyped individuals that are in the pedigree")
    end

    for i in genoID
        push!(ped.setG,i)
    end
    all = Set()
	for i in keys(ped.idMap)
        push!(all,i)
	end
    ped.setNG = setdiff(all,ped.setG)
	j = 1
	for i in ped.setNG
		ped.idMap[i].seqID = j
		j += 1
	end
	numberNonGeno = j - 1
	for i in ped.setG
		ped.idMap[i].seqID = j
		j += 1
	end

	ped.IDs = getIDs(ped) #order of IDs is changed

	return numberNonGeno
end

#for APY
function genoSet!(genoID_file::AbstractString,genoCoreID_file::AbstractString,ped::Pedigree)
    df1 = readtable(genoID_file, eltypes=[String], separator = ' ',header=false)
    for i in df1[:,1]
        push!(ped.setG,i)
	end

    df2 = readtable(genoCoreID_file, eltypes=[String], separator = ' ',header=false)
    for i in df2[:,1]
        push!(ped.setG_core,i)
	end

    all = Set()
	for i in keys(ped.idMap)
        push!(all,i)
	end
    ped.setNG = setdiff(all,ped.setG)
    ped.setG_notcore = setdiff(ped.setG,ped.setG_core)

    j = 1
	for i in ped.setNG
		ped.idMap[i].seqID = j
		j += 1
	end
	numberNonGeno = j - 1

    for i in ped.setG_core
		ped.idMap[i].seqID = j
		j += 1
	end

    for i in ped.setG_notcore
		ped.idMap[i].seqID = j
		j += 1
	end

	numberNonGeno = length(ped.setNG)
	numberGenoCore = length(ped.setG_core)

	ped.IDs=getIDs(ped) #ID are changed

	return (numberNonGeno,numberGenoCore)
end
