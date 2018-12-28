# PedModule

```julia
#Pedigree file
#a 0 0
#b 0 0
#c a b

using PedModule

#Load Pedigree
ped = PedModule.mkPed("ped.txt")

ped.idMap   ##Dict{Any,Any} with 3 entries:
            ##"c" => PedNode(3,"a","b",0.0)
            ##"b" => PedNode(2,"0","0",0.0)
            ##"a" => PedNode(1,"0","0",0.0)

#Calculate the inverse of A (numerator relationship matrix)
Ai = PedModule.AInverse(ped)


#Useful for Single-Step Methods (run it before AInverse(ped))

##Reorder A into 2 groups with the order [others, genotype.ID]
PedModule.genoSet!("genotype.ID",ped)
ped.idMap

##Reorder A into 3 groups with the order
## [others, genotype_core.ID, genotype.ID - genotype_core.ID]
PedModule.genoSet!("genotype.ID","genotype_core.ID",ped)
ped.idMap
```
