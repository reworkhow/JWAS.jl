
![JWAS](assets/JWAS.png)

JWAS is a well-documented software platform based on Julia and an interactive Jupyter notebook for analyses of general
univariate and multivariate Bayesian mixed effects models.  These models are especially useful for, but not limited to,
routine single-trait and multi-trait genomic prediction and genome-wide association studies using either complete or incomplete
genomic data ("single-step" methods). Currently, JWAS provides broad scope of analyses, e.g., a wide collection of Bayesian
methods for whole-genome analyses, including shrinkage estimation and variable selection methods. The features of JWAS include:

* No limitations on fixed effects (e.g. herd-year, age, sex)                                                                    
* Random effects other than markers (e.g. litter, pen)                                  
* Random effects using pedigree information                                                                                
* Random permanent environmental effects  
* Single-trait analyses                                            
* Multi-trait analyses  
* Correlated residuals		
* Correlated random effects
* Correlated marker effects                                                                
* Use of genomic information                                                                                
* Complete genomic data                                      		
* Incomplete genomic data


## Supporting and Citing

We hope the friendly user interface and fast computing speed of JWAS will provide power and convenience for users in both industry
and academia to analyze large datasets. Further, as a well-documented open-source software tool, we hope JWAS will also be used by a
group of active community members, who will contribute to the source code and help maintain the project. Junior scientists can
understand and learn the methodologies for whole-genome analyses by using JWAS and reading the tutorials and source code.

If you would like to help support JWAS, please star the repository on the upper right corner
[here](https://github.com/reworkhow/JWAS.jl) as such statistic will help to demonstrate the
active involvement of the community. If you use JWAS for your research, teaching, or other activities,
we would be grateful if you could cite our work following [this citation guideline](https://github.com/reworkhow/JWAS.jl).

## Get Started:


### Standalone application

A fully self-contained application for JWAS (no installation required) will come out this year.

### Installation

To install julia, please go to the [offical Julia website](https://julialang.org/downloads/).
Please see [platform specific instructions](https://julialang.org/downloads/platform.html)
if you have trouble installing Julia.

To install the package, use the following command inside the Julia REPL (or IJulia Notebook):
```julia
Pkg.add("JWAS")
```

To load the JWAS package, use the following command inside the Julia REPL (or IJulia Notebook):

```julia
using JWAS
```

The command `Pkg.add("JWAS")` will add the registered official `JWAS.jl` and dependencies.

To use the latest/beta features under development, run `Pkg.checkout("JWAS")` to get the
newest unofficial JWAS. Run `Pkg.free("JWAS")` to go back to the offical one.

If you prefer “reproducible research”, an interactive Jupyter notebook interface is available
for Julia (and therefore JWAS). The Jupyter notebook is an open-source web application for creating
and sharing documents that contain live code, equations, visualizations and explanatory text.
To install IJulia, please go to [IJulia](https://github.com/JuliaLang/IJulia.jl).



### the trouble, the error and the new feature

If you have trouble using JWAS, want new featuers or find errors in JWAS, please [open an issue](https://github.com/reworkhow/JWAS.jl/issues) or contact <qtlcheng@ucdavis.edu>.
