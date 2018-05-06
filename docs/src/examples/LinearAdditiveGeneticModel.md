# Bayesian Linear Additive Genetic Model

* link to [Jupyter Notebook](http://nbviewer.jupyter.org/github/reworkhow/
JWAS.jl/blob/master/docs/notebooks_v0.3/2_Linear_Additive_Genetic_Model.ipynb)
for this example

### Step 1: Load Packages

```julia
using JWAS,JWAS.Datasets,CSV,DataFrames
```

### Step 2: Read data

```julia
phenofile  = Datasets.dataset("example","phenotypes.txt")
pedfile    = Datasets.dataset("example","pedigree.txt")

phenotypes = CSV.read(phenofile,delim = ',',header=true)
pedigree   = get_pedigree(pedfile,separator=",",header=true)

head(phenotypes)
```


output:
```julia
6×8 DataFrames.DataFrame
│ Row │ ID │ y1    │ y2   │ y3    │ x1  │ x2 │ x3 │ dam │
├─────┼────┼───────┼──────┼───────┼─────┼────┼────┼─────┤
│ 1   │ a1 │ -0.06 │ 3.58 │ -1.18 │ 0.9 │ 2  │ m  │ 0   │
│ 2   │ a2 │ -0.6  │ 4.9  │ 0.88  │ 0.3 │ 1  │ f  │ 0   │
│ 3   │ a3 │ -2.07 │ 3.19 │ 0.73  │ 0.7 │ 2  │ f  │ 0   │
│ 4   │ a4 │ -2.63 │ 6.97 │ -0.83 │ 0.6 │ 1  │ m  │ a2  │
│ 5   │ a5 │ 2.31  │ 3.5  │ -1.52 │ 0.4 │ 2  │ m  │ a2  │
│ 6   │ a6 │ 0.93  │ 4.87 │ -0.01 │ 5.0 │ 2  │ f  │ a3  │
```


## Univariate Linear Additive Genetic Model

### Step 3: Build Model Equations

```julia
model_equation = "y1 = intercept + x1*x3 + x2 + x3 + ID + dam"
R     = 1.0
model = build_model(model_equation,R);
```

### Step 4: Set Factors or Covariate
```julia
set_covariate("x1");
```

### Step 5: Set Random or Fixed Effects
```julia
G1 = 1.0
G2 = eye(2)
set_random(model,"x2",G1)
set_random(model,"ID dam",pedigree,G2)
```

### Step 6: Run Bayesian Analysis
```julia
out=runMCMC(model,phenotypes)
```


## Multivariate Linear Additive Genetic Model

### Step 3: Build Model Equations

```julia
model_equation = "y1 = intercept + x1 + x3 + ID + dam
                  y2 = intercept + x1 + x2 + x3 + ID  
                  y3 = intercept + x1 + x1*x3 + x2 + ID"

R     = eye(3)
model = build_model(model_equation,R);
```

### Step 4: Set Factors or Covariate
```julia
set_covariate("x1");
```

### Step 5: Set Random or Fixed Effects
```julia
G1 = eye(2)
G2 = eye(4)
set_random(model,"x2",G1)
set_random(model,"ID dam",pedigree,G2)
```

### Step 6: Run Bayesian Analysis
```julia
out=runMCMC(model,phenotypes)
```
