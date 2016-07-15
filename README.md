JWAS.jl
=========

[![Build Status](https://travis-ci.org/reworkhow/JWAS.jl.svg?branch=master)](https://travis-ci.org/reworkhow/JWAS.jl)

JWAS.jl is an open-source software tool written in Julia for Bayesian multiple regression methods applied to genome-wide association studies and genomic prediction.

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.clone("https://github.com/reworkhow/JWAS.jl.git")`
* ~~**Documentation**: [available here](http://jwasjl.readthedocs.org/en/latest/)~~
* **Examples**: [available here](http://nbviewer.jupyter.org/github/reworkhow/JWAS.jl/tree/master/test/)

<figure class="highlight"><pre><code class="language-shell" data-lang="shell">
<mark style="background-color:red;"><big>JWAS.jl</big></mark>

├──────── <mark style="background-color:orange;">PedModule.jl</mark>

├──────── <mark style="background-color:orange;">ST.jl</mark>
           ├── <i>build_model</i>
           ├── <i>set_covariate</i>
           ├── <i>set_random</i>
           ├── <i>get_pedigree</i>
           ├── <i>add_markers</i>
           ├── <i>outputMCMCsamples</i>
           ├── <i>showMME</i>
           ├── <i>solve</i>
           └── <i>runMCMC</i>

├──────── <mark style="background-color:orange;">MT.jl</mark>
           ├── <i>MT.build_model</i>
           ├── <i>MT.set_covariate</i>
           ├── <i>MT.set_random</i>
           ├── <i>MT.get_pedigree</i>
           ├── <i>MT.add_markers</i>
           ├── <i>MT.showMME</i>
           ├── <i>MT.solve</i>
           └── <i>MT.runMCMC</i>

├──────── <mark style="background-color:orange;">???.jl</mark>
</code>
</pre>
</figure>
