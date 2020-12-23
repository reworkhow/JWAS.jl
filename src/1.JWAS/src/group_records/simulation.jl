using StatsBase,Statistics, Distributions,Printf,Random,DelimitedFiles,InteractiveUtils,DataFrames,CSV,SparseArrays,LinearAlgebra
# pedigree information: Ind: individual number, Sire: individual's father, Dam: individual's mother, Sex: 1=male, 0=female, generation: 1 to 8, Breeding: 1=selected as breeding animal, 0 = otherwise
ped = DataFrame(Ind = Int[], Sire = Int[], Dam = Int[], Sex = Int[], generation = Int[], Breeding = Int[])
# Genration 1
for i in 1:30  # 30 sires
    push!(ped,(size(ped)[1]+1, 0, 0, 1, 1, 1))
end

for i in 1:20  # 20 dams
    push!(ped,(size(ped)[1]+1, 0, 0, 0, 1, 1))
end

# Generation 2 to 8
Random.seed!(101)
for g in 2:8
    #select breeding animals from last generation
    sire_select = sample(findall(ped.Sex[ped.generation .== (g-1)] .== 1), 30, replace = false) .+ (50*(g>2) + (g-3)*3600*(g>2))
    dam_select = sample(findall(ped.Sex[ped.generation .== (g-1)] .== 0), 20, replace = false) .+ (50*(g>2) + (g-3)*3600*(g>2))
    ped.Breeding[sire_select] .= 1
    ped.Breeding[dam_select] .= 1

    # individual information
    for i in sire_select
        for j in dam_select
            for k in 1:6
                push!(ped, (size(ped)[1]+1, i, j, rand((0,1)), g, 0))
            end
        end
    end
end

litter = zeros(Int, 50)
for i in 1:Int((nrow(ped)-50)/6)
    append!(litter, ones(Int, 6)*i)
end
ped.litter = litter

function mkmat_incidence_factor(yID,uID)
    Z = spzeros(length(yID),length(uID))

    uIDdict = Dict()
    for (index,id) in enumerate(uID)
        uIDdict[id]=index
    end

    rowi = 1
    for id in yID
        if haskey(uIDdict,id)
            index = uIDdict[id]
        else
            error(id, " is not found!")
        end
        Z[rowi,index]=1
        rowi = rowi+1
    end
    return Z
end

## simulate covariates and phenotype
Random.seed!(102)
# only last five generations used for data analysis
ped_ana = ped[ped.generation .> 3, :]

# litter effect, variance = 40
l_eff = randn(length(unique(ped_ana.litter))) .* (40^0.5)
l_mat = mkmat_incidence_factor(ped_ana.litter, unique(ped_ana.litter))
l_val = l_mat * l_eff

# pen allocation and pen effect, variance = 40
pen_size = 12
pen_n = Int(nrow(ped_ana)/pen_size)
# 1) two litters were assigned to the same pen
pen_allocation = repeat(append!(shuffle(1:pen_n),shuffle(1:pen_n)), inner = 6)
ped_ana.pen1 = pen_allocation
# 2) a litter was randomly divided into two sub-litters of size 3, which were randomly allocated to pens
pen_all = []
for i in 1:4
    append!(pen_all, shuffle(1:pen_n))
end
pen_all = repeat(pen_all, inner = 3)
pen_allocation = []
for i in 1:Int(nrow(ped_ana)/6)
    append!(pen_allocation, shuffle(pen_all[((i-1)*6 + 1):(i*6)]))
end
ped_ana.pen2 = pen_allocation

# 3) complete randomization
pen_allocation = shuffle(repeat(1:pen_n, inner = 12))
ped_ana.pen3 = pen_allocation

pen_eff = randn(pen_n) .* (40^0.5)
pen_mat = mkmat_incidence_factor(ped_ana.pen1, unique(ped_ana.pen1))
pen_val = pen_mat * pen_eff

# Additive genetic effects ?

# Residual, variance = 200 and 280
res1 = randn(nrow(ped_ana)) .* (200^0.5)
res2 = randn(nrow(ped_ana)) .* (280^0.5)

# Phenotypic values
pheno1 = l_val + pen_val + gen_eff + res1
pheno2 = gen_eff + res2

# keep 80% of the observations (but keep all the breeding animals)
deleterows!(ped_ana, sort(sample(findall(ped_ana.Breeding .== 0),Int(nrow(ped_ana)*0.2), replace = false)))
