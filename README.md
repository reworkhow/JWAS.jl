# JWAS.jl

[![Build Status](https://travis-ci.com/reworkhow/JWAS.jl.svg?branch=master)](https://travis-ci.org/reworkhow/JWAS.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://reworkhow.github.io/JWAS.jl/latest)
<!---[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://reworkhow.github.io/JWAS.jl/stable)--->

![JWAS](docs/src/assets/JWAS.png)

JWAS.jl is an open-source software tool written in Julia for Bayesian multiple regression methods applied to genomic prediction and genome-wide association studies.

* **Homepage**: [QTL.rocks](https://QTL.rocks)
* **Discussion group**: [available here](https://groups.io/g/qtlrocks)
* **Installation**: at the Julia REPL, `using Pkg; Pkg.add("JWAS")`
* **Documentation**: [available here](http://reworkhow.github.io/JWAS.jl/latest/)
* **Examples**: [available here](https://github.com/reworkhow/JWAS.jl/wiki)
* **Authors**: [Hao Cheng](http://qtl.rocks/chenglab.html), [Rohan Fernando](http://www.ans.iastate.edu/people/rohan-fernando), [Dorian Garrick](http://www.massey.ac.nz/massey/expertise/profile.cfm?stref=600130)



### Structure of JWAS

```
JWAS.jl

├────── get_pedigree
├────── get_genotypes
├────── build_model
├────── set_covariate
├────── set_random
├────── solve
├────── runMCMC
└────── ...
```

### Help

1. Show this README file in REPL or notebook using `?JWAS`
2. For help on a specific function above, type ? followed by its name, e.g. `?runMCMC` and press enter.
3. Run `Pkg.add(PackageSpec(name="JWAS", rev="master"))` to get the newest unofficial JWAS. Run `Pkg.free("JWAS")` to go back to the official one.

### Examples

 * [Single Trait Analysis](https://github.com/reworkhow/JWAS.jl/wiki/Single-Trait-Analysis)

 * [Multiple Trait Analysis](https://github.com/reworkhow/JWAS.jl/wiki/Multiple-Trait-Analysis)

 * [Repeated Measures](https://github.com/reworkhow/JWAS.jl/wiki/Repeated-Measures)

 * [Single Step Analysis](https://github.com/reworkhow/JWAS.jl/wiki/Single-Step-Analysis)

 * [Categorical Trait Analysis](https://github.com/reworkhow/JWAS.jl/wiki/Categorical-Trait-Analysis)

 * [Censored Trait Analysis](https://github.com/reworkhow/JWAS.jl/wiki/Censored-Trait-Analysis)

 * [Multi-class Bayesian Analysis](https://github.com/reworkhow/JWAS.jl/wiki/Multiclass-Bayesian-analysis)

 * [Neural Networks](https://github.com/reworkhow/JWAS.jl/wiki/Neural-Networks)

 * [Cross Validation](https://github.com/reworkhow/JWAS.jl/wiki/Cross-Validation)	

 * [Genome Wide Association Study](https://github.com/reworkhow/JWAS.jl/wiki/GWAS)	

 * [Integrating Phenotypic Causal Networks in GWAS](https://github.com/reworkhow/JWAS.jl/wiki/Integrating-Phenotypic-Causal-Networks-in-GWAS)

 * [single trait and multiple trait GBLUP by providing the relationship matrix directly](https://github.com/reworkhow/JWAS.jl/wiki/single-trait-and-multiple-trait-GBLUP-by-providing-the-relationship-matrix-directly)

 * [Description of Mixed Effects Model](https://github.com/reworkhow/JWAS.jl/wiki/Description-of-Mixed-Effects-Model)

 * [Quality Control of Genotypes](https://github.com/reworkhow/JWAS.jl/wiki/Quality-Control-of-Genotypes)
