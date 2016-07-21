# PedModule

[![Build Status](https://travis-ci.org/reworkhow/PedModule.jl.svg?branch=master)](https://travis-ci.org/reworkhow/PedModule.jl)

###pedigree file

```bash
a 0 0
b 0 0
c a b
```

####Quick-start

```Julia
using PedModule

#input pedigree
ped = PedModule.mkPed("ped.txt")

ped.idMap   ##Dict{Any,Any} with 3 entries:
            ##"c" => PedNode(3,"a","b",0.0)
            ##"b" => PedNode(2,"0","0",0.0)
            ##"a" => PedNode(1,"0","0",0.0)
            
#calculate A inverse
Ai = PedModule.AInverse(ped)

##reorder A into 2 groups (specifically for single step methods) in the order [others,genotype.ID]
##run it before PedModule.AInverse(ped)
PedModule.genoSet!("genotype.ID",ped)
ped.idMap

##reorder A into 3 groups (specifically for single step methods) in the order [others,genotype_core.ID,genotype.ID-genotype_core.ID]
##run it before PedModule.AInverse(ped)
PedModule.genoSet!("genotype.ID","genotype_core.ID",ped)
ped.idMap
```

####More

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.clone("https://github.com/QTL-rocks/PedModule.jl.git")`
* **Authors**: [Rohan Fernando](http://www.ans.iastate.edu/faculty/index.php?id=rohan), [Hao Cheng](http://reworkhow.github.io)
