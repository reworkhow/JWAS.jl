# Workflow

## Check Julia version

## Read data

```julia
data = CSV.read("data.txt")
```

## Build Model Equations

```julia
model=build_model("y=x1+x2+x3+x4")
```

- link to [`build_model`](@ref)

## Set Factors or Covariate
```julia
model=build_model("y=x1+x2+x3+x4")
```

- link to [`build_model`](@ref)


## Set Random or Fixed Effects
```julia
model=build_model("y=x1+x2+x3+x4")
```

- link to [`build_model`](@ref)


## Use Pedigree Information
```julia
model=build_model("y=x1+x2+x3+x4")
```

- link to [`build_model`](@ref)


## Use Genomic Information
```julia
model=build_model("y=x1+x2+x3+x4")
```

- link to [`build_model`](@ref)


- link to [Workflow](@ref)
