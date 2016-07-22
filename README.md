# JWAS.jl

[![Build Status](https://travis-ci.org/reworkhow/JWAS.jl.svg?branch=master)](https://travis-ci.org/reworkhow/JWAS.jl)

JWAS.jl is an open-source software tool written in Julia for Bayesian multiple regression methods applied to genome-wide association studies and genomic prediction.

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.add("JWAS")`
* ~~**Documentation**: [available here](http://jwasjl.readthedocs.org/en/latest/)~~
* **Examples**: [available here](http://nbviewer.jupyter.org/github/reworkhow/JWAS.jl/tree/master/test/)



### Structure of JWAS Module


```
JWAS.jl

├──────── PedModule.jl

├──────── ST.jl
           ├── build_model
           ├── set_covariate
           ├── set_random
           ├── get_pedigree
           ├── add_markers
           ├── outputMCMCsamples
           ├── showMME
           ├── solve
           └── runMCMC

├──────── MT.jl
           ├── MT.build_model
           ├── MT.set_covariate
           ├── MT.set_random
           ├── MT.get_pedigree
           ├── MT.add_markers
           ├── MT.showMME
           ├── MT.solve
           └── MT.runMCMC

├──────── QTL.jl
           ├── get_additive_genetic_variances
           └── get_breeding_values
           
├──────── Datasets.jl

```

### Help

1. Show this README file in REPL or notebook using `?JWAS`
2. For help on a specific function above, type ? followed by its name, e.g. `?runMCMC` and press enter.
