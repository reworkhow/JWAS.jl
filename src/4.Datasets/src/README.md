# Dataset Simulation Scripts

This folder contains standalone data-generation scripts used to create example
datasets for JWAS tutorials and local experimentation. These scripts are not
part of the package runtime API. They are entry points with side effects: they
read or write files in the current working directory.

## Files

### `simulation.jl`

Self-contained XSim-based generator for a synthetic multi-trait example.

What it does:
- builds a small genome with XSim
- simulates founders and several pedigree cohorts
- exports `pedigree.txt` and `genotypes.txt`
- simulates three-trait phenotypes with a mix of breeding values and
  structured non-genetic effects
- exports `phenotypes.txt`

Generated outputs:
- `output.txt.ped`
- `output.txt.gen`
- `pedigree.txt`
- `genotypes.txt`
- `phenotypes.txt`

Notes:
- only the first 1000 SNPs are written to `genotypes.txt`
- the remaining loci act as extra background/polygenic signal in the phenotype simulation

### `simulation_realgeno.jl`

Phenotype simulator that starts from an existing genotype matrix and pedigree.

Required input files in the current working directory:
- `genotypes.csv`
- `pedigree.csv`

What it does:
- samples QTL from the marker panel
- simulates correlated breeding values for three traits
- adds a pedigree-based maternal effect to the first trait
- adds fixed and random non-genetic effects
- writes `phenotypes.csv`

Generated output:
- `phenotypes.csv`

## Usage

Run these scripts from a directory where you want the output files to be created
and where any required input files are already present.

Examples:

```bash
julia --project=. src/4.Datasets/src/simulation.jl
julia --project=. src/4.Datasets/src/simulation_realgeno.jl
```

Because these scripts write files directly, avoid running them from the
repository root unless you intend to overwrite local scratch data there.
