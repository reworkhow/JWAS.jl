# Workflow

A step by step workflow for how to run JWAS is shown in this section. Given the data and model equations, different types of models
as shown in the table below can be  used. The "X" denotes the type of available data. The "A <= B" denotes that A individuals is a subset of B individuals.  


| Linear Mixed Models (LMM) |    phenotypes | pedigree  | genotypes | notes                  |
| ------------------------- |:-------------:|:---------:|:---------:|:----------------------:|
| Conventional LMM          |              X|           |           |                        |
| Pedigree-based LMM        |              X|         X |           | phenotypes <= pedigree |
| Complete Genomic LMM      |              X|     maybe |          X| phenotypes <= genotypes|
| Incomplete Genomic LMM    |              X|         X |          X| phenotypes <= pedigree and genotypes <= pedigree|


* **Incomplete Genoimc LMM** is also called "single-step" methods in animal breeding.

* Pedigree information may be used in **Complete Genomic LMM** as a seperate polygenic term to account for genetic variance not explained by the genomic information (e.g., SNPs).

* Note that **Pedigree-based LMM** and **Complete Genomic LMM** are special cases of **Incomplete Genomic LMM**.


## Get Data Ready

By default, input data files are comma-separated values (CSV) files, where each line of the file is a data record, and each record consists of one or more fields, separated by **commas**. Other field separators such as space (' ') or tab ('\t') can also be used if you supply the keyword argument, e.g, **CSV.read(...,delim='\t')** or **add_genotypes(...,separator='\t')**

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

## Load packages
```julia
using JWAS, CSV, DataFrames
```


## Read data

```julia
data = CSV.read("data.txt")
head(data)
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
## Build Model Equations

```julia
model_equation = "y1 = x1 + x3;
                  y2 = x1 + x2 +x3;  
                  y2 = x1 + x1*x3 + x2"
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
set_random("x2",0.1)
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
