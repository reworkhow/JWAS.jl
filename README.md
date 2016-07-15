# JWAS.jl

[![Build Status](https://travis-ci.org/reworkhow/JWAS.jl.svg?branch=master)](https://travis-ci.org/reworkhow/JWAS.jl)

JWAS.jl is an open-source software tool written in Julia for Bayesian multiple regression methods applied to genome-wide association studies and genomic prediction.

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.clone("https://github.com/reworkhow/JWAS.jl.git")`
* ~~**Documentation**: [available here](http://jwasjl.readthedocs.org/en/latest/)~~
* **Examples**: [available here](http://nbviewer.jupyter.org/github/reworkhow/JWAS.jl/tree/master/test/)


### Structure of JWAS Module

```
JWAS.j

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
```

### Troubleshooting

1. Get help about functions above through **?foo** or **@doc(foo)**

