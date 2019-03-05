# JWAS.jl

[![Build Status](https://travis-ci.org/reworkhow/JWAS.jl.svg?branch=master)](https://travis-ci.org/reworkhow/JWAS.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://reworkhow.github.io/JWAS.jl/latest)
<!---[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://reworkhow.github.io/JWAS.jl/stable)--->


JWAS.jl is an open-source software tool written in Julia for Bayesian multiple regression methods applied to genome-wide association studies and genomic prediction.

* **Homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `using Pkg; Pkg.add("JWAS")`
* **Documentation**: [available here](http://reworkhow.github.io/JWAS.jl/latest/)
* **Notebooks**: [available here](http://nbviewer.jupyter.org/github/reworkhow/JWAS.jl/tree/master/docs/notebooks/0_index.ipynb)
* **Authors**: [Hao Cheng](http://qtl.rocks/chenglab.html), [Rohan Fernando](http://www.ans.iastate.edu/people/rohan-fernando), [Dorian Garrick](http://www.massey.ac.nz/massey/expertise/profile.cfm?stref=600130)



### Structure of JWAS Module



```
JWAS.jl

├────── build_model
├────── set_covariate
├────── set_random
├────── get_pedigree
├────── add_genotypes
├────── outputMCMCsamples
├────── showMME
├────── solve
└────── runMCMC

├──────── PedModule.jl

├──────── Datasets.jl

├──────── misc.jl
           ├── get_breeding_values
           ├── get_additive_genetic_variances
           ├── get_heritability
           ├── get_correlations
           └── report

```

### Help

1. Show this README file in REPL or notebook using `?JWAS`
2. For help on a specific function above, type ? followed by its name, e.g. `?runMCMC` and press enter.
3. Run `Pkg.add(PackageSpec(name="JWAS", rev="master"))` to get the newest unofficial JWAS. Run `Pkg.free("JWAS")` to go back to the offical one.
