# Workflow

A step by step workflow for how to run JWAS is shown in this section. Given the data and model equations, different types of models
as shown in the table below can be  used. The "X" denotes the type of available data. The "A <= B" denotes that A individuals is a subset of B individuals.  


| Linear Mixed Models (LMM) |    phenotypes | pedigree  | genotypes | notes                  |
| ------------------------- |:-------------:|:---------:|:---------:|:----------------------:|
| conventional LMM          |              X|           |           |                        |
| pedigree-based LMM        |              X|         X |           | phenotypes <= pedigree |
| complete genomic LMM      |              X|     maybe |          X| phenotypes <= genotypes|
| incomplete genomic LMM    |              X|         X |          X| phenotypes <= pedigree and genotypes <= pedigree|


* Pedigree information may be fitted in complete genomic LMM as a seperate polygenic term to account for genetic variance not
explained by the genomic data. pedigree-based LMM  and complete genomic LMM are sepcial cases of incomplete genomic LMM.

* Note that pedigree-based LMM  and complete genomic LMM are sepcial cases of incomplete genomic LMM.


## Get Data Ready

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

<p>Click on the buttons inside the tabbed menu to see the data format:</p>

<div class="tab">
  <button class="tablinks" onclick="openCity(event, 'phenotypes')">phenotypes.txt</button>
  <button class="tablinks" onclick="openCity(event, 'pedigree')">pedigree.txt</button>
  <button class="tablinks" onclick="openCity(event, 'genotypes')">genotypes.txt</button>
</div>

<div id="phenotypes" class="tabcontent">
  <p>id,y1,y2,y3,x1,x2,x3,x4,dam<br />
     a1,y1,y2,y3,x1,x2,x3,x4,dam<br />
     a2,y1,y2,y3,x1,x2,x3,x4,dam<br />
     a3,y1,y2,y3,x1,x2,x3,x4,dam<br />
     a4,y1,y2,y3,x1,x2,x3,x4,dam<br />
     a5,y1,y2,y3,x1,x2,x3,x4,dam<br />
     a6,y1,y2,y3,x1,x2,x3,x4,dam<br />
     a7,y1,y2,y3,x1,x2,x3,x4,dam<br />
     a10,y1,y2,y3,x1,x2,x3,x4,dam<br />
     a,y1,y2,y3,x1,x2,x3,x4,dam<br />

  </p>
</div>

<div id="pedigree" class="tabcontent">
<p>id,sire,dam</p>
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
<p>id,m1,m2,m3,m4,m5</p>
<p>a1,1,2,1,1,0</p>
<p>a2,2,1,1,1,1</p>
<p>a3,1,1,0,1,1</p>
<p>a4,0,1,2,1,0</p>
<p>a5,0,1,1,1,1</p>
<p>a6,1,2,0,0,0</p>
<p>a7,2,1,1,1,1</p>
<p>a8,1,2,1,0,1</p>
<p>a9,1,2,2,0,0</p>
<p>a10,1,2,1,1,1</p>
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
