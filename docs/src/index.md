
![JWAS](assets/JWAS.png)

JWAS is a well-documented software platform based on Julia and an interactive Jupyter notebook for analyses of general
univariate and multivariate Bayesian mixed effects models.  These models are especially useful for, but not limited to,
routine single-trait and multi-trait genomic prediction and genome-wide association studies using either complete or incomplete
genomic data ("single-step" methods). Currently, JWAS provides broad scope of analyses, e.g., a wide collection of Bayesian
methods for whole-genome analyses, including shrinkage estimation and variable selection methods. The features of JWAS include:

* Univariate (single-trait) analysis
* Multivariate (multi-trait) analysis  
* No limitations on fixed effects (e.g., herd, year, age, sex)
* Random effects other than markers (e.g., litter, pen)                                  
* Random effects using pedigree information
  - Additive genetic effects
  - Maternal effects
* Random permanent environmental effects  
* Correlated residuals		
* Correlated random effects
* Unknown (or known) variance components
* Use of genomic information
  - Complete genomic data                                      		
  - Incomplete genomic data (singe-step)


## Supporting and Citing

We hope the friendly user interface and fast computing speed of JWAS will provide power and convenience for users in both industry
and academia to analyze large datasets. Further, as a well-documented open-source software tool, we hope JWAS will also be used by a
group of active community members, who will contribute to the source code and help maintain the project. Junior scientists can
understand and learn the methodologies for whole-genome analyses by using JWAS and reading the tutorials and source code.

If you would like to help support JWAS, please star the repository on the upper right corner
[here](https://github.com/reworkhow/JWAS.jl) as such statistic will help to demonstrate the
active involvement of the community. If you use JWAS for your research, teaching, or other activities,
we would be grateful if you could cite our work following [Citing](@ref).


## The trouble, the error and the new feature

If you have trouble using JWAS, want new features or find errors in JWAS, please [open an issue](https://github.com/reworkhow/JWAS.jl/issues) or contact <qtlcheng@ucdavis.edu>.

## Tutorials

### Theory
```@contents
Pages = [
  "theory/theory.md"
]
Depth = 3
```

### Manual
```@contents
Pages = [
  "manual/getstarted.md",
  "manual/workflow.md",
  "manual/public.md",
  "manual/internals.md",
]
Depth = 3
```

### Examples
```@contents
Pages = [
  "examples/examples.md"
]
Depth = 2
```
