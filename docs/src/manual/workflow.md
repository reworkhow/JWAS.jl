# Workflow

A step by step workflow for how to run JWAS is shown in this section. The workflow below
is used to demonstrate a three-trait Bayesian linear mixed model fitting fixed effects (x1, x3),
random effects (x2), direct genetic effects (ID), maternal genetic effects (dam) and
genomic information.

```@contents
Pages = [
  "workflow.md"
]
Depth = 2
```


## Available Models

Given the data and model equations, several different types of models are available in JWAS as shown below. In the table below, "X"
denotes the type of available data, and "Y<=A" denotes that Y individuals is a subset of A individuals.  


| Linear Mixed Models (LMM) | phenotypes (Y)| pedigree (A)| genotypes (G)| notes     |
| :-----------------------: |:-------------:|:-----------:|:------------:|:---------:|
| Conventional LMM          |              X|             |              |           |
| Pedigree-based LMM        |              X|         X   |              | Y<=A      |
| Complete Genomic LMM      |              X|     optional|             X| Y<=G      |
| Incomplete Genomic LMM    |              X|         X   |             X| Y<=A,G<=A |


!!! note

    - **Incomplete Genomic LMM** is also called "single-step" methods in animal breeding.

    - Pedigree information may be used in **Complete Genomic LMM** for extra polygenic effects to account for genetic variance not explained by the genomic data (e.g., SNPs).

    - **Pedigree-based LMM** (none of the individuals in the pedigree are genotyped) and **Complete Genomic LMM** (all individuals in the pedigree are genotyped) are
      special cases of **Incomplete Genomic LMM** (part of the individuals in the pedigree are genotyped).


## Get Data Ready

By default, input data files are comma-separated values (CSV) files, where each line of the file consists of one or more fields, separated by **commas**. Other field separators such as space (' ') or tab ('\t') can be used if you supply the keyword argument, e.g, `CSV.read(...,delim='\t')` or `add_genotypes(...,separator='\t')`

Click on the buttons inside the tabbed menu to see the data:


```@raw html
<head>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
body {font-family: Arial;}

/* Style the tab */
.tab {
    overflow: hidden;
    border: 1px solid #ccc;
    background-color: #f1f1f1;
}

/* Style the buttons inside the tab */
.tab button {
    background-color: inherit;
    float: left;
    border: none;
    outline: none;
    cursor: pointer;
    padding: 14px 16px;
    transition: 0.3s;
    font-size: 17px;
}

/* Change background color of buttons on hover */
.tab button:hover {
    background-color: #ddd;
}

/* Create an active/current tablink class */
.tab button.active {
    background-color: #ccc;
}

/* Style the tab content */
.tabcontent {
    display: none;
    padding: 6px 12px;
    border: 1px solid #ccc;
    border-top: none;
}
</style>
</head>
<body>

<div class="tab">
  <button class="tablinks" onclick="openCity(event, 'phenotypes')">phenotypes.txt</button>
  <button class="tablinks" onclick="openCity(event, 'pedigree')">pedigree.txt</button>
  <button class="tablinks" onclick="openCity(event, 'genotypes')">genotypes.txt</button>
  <button class="tablinks" onclick="openCity(event, 'map file')">map.txt</button>
</div>

<div id="phenotypes" class="tabcontent">
<p>ID,y1,y2,y3,x1,x2,x3,dam</p>
<p>a1,-0.06,3.58,-1.18,0.9,2,m,0</p>
<p>a2,-0.6,4.9,0.88,0.3,1,f,0</p>
<p>a3,-2.07,3.19,0.73,0.7,2,f,0</p>
<p>a4,-2.63,6.97,-0.83,0.6,1,m,a2</p>
<p>a5,2.31,3.5,-1.52,0.4,2,m,a2</p>
<p>a6,0.93,4.87,-0.01,05,2,f,a3</p>
<p>a7,-0.69,3.1,-1.47,0.5,2,f,a3</p>
<p>a8,-4.69,7.31,-1.09,0.3,2,m,a6</p>
<p>a9,-2.81,7.18,0.76,0.4,2,m,a6</p>
<p>a10,1.92,1.78,-0.88,0.2,1,m,a7</p>
</div>

<div id="pedigree" class="tabcontent">
<p>ID,Sire,Dam</p>
<p>a1,0,0</p>
<p>a2,0,0</p>
<p>a3,0,0</p>
<p>a4,a1,a2</p>
<p>a5,a1,a2</p>
<p>a6,a1,a3</p>
<p>a7,a1,a3</p>
<p>a8,a4,a6</p>
<p>a9,a4,a6</p>
<p>a10,a5,a7</p>
</div>

<div id="genotypes" class="tabcontent">
<p>ID,m1,m2,m3,m4,m5</p>
<p>a1,1,2,1,1,0</p>
<p>a2,2,1,1,1,1</p>
<p>a3,1,1,0,1,1</p>
<p>a4,2,2,0,1,0</p>
<p>a5,1,1,2,1,1</p>
<p>a6,2,1,0,0,0</p>
<p>a7,0,2,1,2,1</p>
<p>a8,2,2,0,0,0</p>
<p>a9,2,1,0,1,0</p>
<p>a10,0,2,2,2,1</p>
</div>

<div id="map file" class="tabcontent">
<p>markerID,chromosome,position</p>
<p>m1,1,16977</p>
<p>m2,1,434311</p>
<p>m3,1,1025513</p>
<p>m4,2,70350</p>
<p>m5,2,101135</p>
</div>


<script>
function openCity(evt, cityName) {
    var i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
    }
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
    }
    document.getElementById(cityName).style.display = "block";
    evt.currentTarget.className += " active";
}
</script>
</body>
```

## Step 1: Load Packages

```julia
using JWAS,CSV,DataFrames
```
---
The `JWAS` package is loaded, as well as the `CSV` and `DataFrame` packages for reading text files.


## Step 2: Read data

```julia
phenotypes = CSV.read("phenotypes.txt",delim = ',',header=true)
pedigree   = get_pedigree("pedigree.txt",separator=",",header=true)
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

- link to documentation for [`get_pedigree`](@ref)

---
The **phenotypic data** is read on line 1, and the **pedigree data** is read on line 2. On line 3, the first several rows of data are shown.

## Step 3: Build Model Equations

```julia
model_equation = "y1 = intercept + x1 + x3 + ID + dam
                  y2 = intercept + x1 + x2 + x3 + ID  
                  y3 = intercept + x1 + x1*x3 + x2 + ID"
model=build_model(model_equation, R)
```

- link to documentation for [`build_model`](@ref)

---
The non-genomic part of the model equation for a 3-trait analysis is defined on the first 3 lines.
* The effects fitted in the model for trait `y1` are the intercept, `x1`, `x3`, direct genetic effects (`ID`) and maternal genetic effects (`dam`).
* The effects fitted in the model for trait `y2` are the intercept, `x1`, `x2`, `x3` and direct genetic effects (`ID`).
* The effects fitted in the model for trait `y3` are the intercept, `x1`, the interaction between `x1` and `x3`, `x2` and direct genetic effects (`ID`).

On the last line, the model is built given the model equation and residual variance `R` (a 3x3 matrix).
By default, all effects are treated as fixed and classed as factors (categorical variables)
rather than covariates (quantitative variables).

## Step 4: Set Factors or Covariate
```julia
set_covariate(model,"x1")
```

- link to documentation for [`set_covariate`](@ref)

---
On line 1, the effect `x1` is defined to be a covariate rather than class effect.

## Step 5: Set Random or Fixed Effects
```julia
set_random(model,"x2",G1)
set_random(model,"ID dam",pedigree,G2)
```
- link to documentation for [`set_random`](@ref)

---
On line 1, the `x2` class effect is defined as random with variance `G1`(a 2x2 matrix). On line 2, direct genetic effects and
maternal genetic effects are fitted as `ID` and `dam` with `G2` (a 4x4 matrix) and the inverse of the numerator relationship matrix defined from pedigree.

## Step 6: Use Genomic Information
```julia
add_genotypes(model,"genotypes.txt",G3,separator=',')
```

- link to documentation for [`add_genotypes`](@ref)

---
On line 1, the genomic part of the model is defined with the genotype file and variance `G3` (a 3x3 matrix).

## Step 7: Run Bayesian Analysis
```julia
outputMCMCsamples(model,"x2")
out=runMCMC(model,phenotypes,methods="BayesC",output_samples_frequency=100)
```

- link to documentation for [`outputMCMCsamples`](@ref)
- link to documentation for [`runMCMC`](@ref)

---
On line 1, MCMC samples from `runMCMC` for `x2` is saved to a file, where each row represents one sample from the MCMC.
On line 2, a multi-trait BayesC analysis is performed with `model` and `phenotypes` as had been defined in step 1-6.
MCMC samples for marker effects, location parameters specified on line 1, and all variance components from this analysis
are saved every `output_samples_frequency` iterations to files.

---
Several steps above can be skipped if no related information is available, e.g., step 6 is skipped
for pedigree-based LMM. Several detailed examples are available in the examples section. Here is the link
to documentation for all [Public functions](@ref).

## check results

Posterior means of location parameters, most variance components, and marker effects are saved in `out`.
They can be listed and obtained as
```julia
keys(out)

# output:
#
# Base.KeyIterator for a Dict{Any,Any} with 7 entries. Keys:
#   "Posterior mean of polygenic effects covariance matrix"
#   "Model frequency"
#   "Posterior mean of residual covariance matrix"
#   "Posterior mean of marker effects"
#   "Posterior mean of marker effects covariance matrix"
#   "Posterior mean of location parameters"
#   "Posterior mean of Pi"

out["Posterior mean of residual covariance matrix"]

# output:
#
# 3×3 Array{Float64,2}:
#   0.674651   -0.103877   0.0834044
#  -0.103877    0.828135  -0.121798
#   0.0834044  -0.121798   0.720751

```

MCMC samples for marker effects, location parameters specified in step 7, and all variance components are saved to text
files in your working directory. They can be obtained as

```julia
res=readdlm("MCMC_samples_marker_effects_y1.txt",',',header=true)
```
