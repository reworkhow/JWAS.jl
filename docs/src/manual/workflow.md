# Workflow



## Data format

```
#data.txt
id,y1,y2,y3,x1,x2,x3,x4,dam

#pedigree
id,sire,dam

#genotype
id,m1,m2,m3,m4,m5

```

## Read data

```julia
data = CSV.read("data.txt")
```


## Build Model Equations

```julia
model_equation = "y1 = x1 + x2 + x3 + x4;
                  y2 = x1 + x2 + x3*x4"
model=build_model(model_equation)
```

- link to [`build_model`](@ref)

## Set Factors or Covariate
```julia
set_covariate("x1")
```

- link to [`set_covariate`](@ref)


## Set Random or Fixed Effects
```julia
set_random("x2",)
```

- link to [`set_random`](@ref)


## Use Pedigree Information
```julia
ped=get_pedigree("pedigree.txt")
```

- link to [`get_pedigree`](@ref)


## Use Genomic Information
```julia
add_genotypes(model,"genotypes.txt")
```

- link to [`add_genotypes`](@ref)


- link to [Workflow](@ref)
