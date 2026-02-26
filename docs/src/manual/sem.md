# SEM: Beginner to Advanced

This page is a single, layered guide to structural equation model (SEM) usage in JWAS.
It is organized from practical setup to exact reference behavior to implementation-level notes.

```@contents
Pages = ["sem.md"]
Depth = 2
```

## Quick Index

- [Beginner Quick Start](#Beginner-Quick-Start)
- [Reference](#Reference)
- [Deep Technical Notes](#Deep-Technical-Notes)
- [Troubleshooting](#Troubleshooting)

## Beginner Quick Start

Use SEM in JWAS when you have multiple phenotypes and a directed causal order among traits.
JWAS samples structural coefficients during MCMC and reports direct, indirect, and overall marker effects.

### Requirements and Constraints

- SEM is available in multi-trait analysis only.
- `causal_structure` must be lower triangular.
- Matrix interpretation: column index affects row index.
- When SEM is enabled, JWAS enforces residual covariance constraint and disables missing phenotype imputation for identifiability.

### End-to-End Example

```julia
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets

# Step 1: Read demo data
phenofile  = dataset("phenotypes.txt", dataset_name="demo_7animals")
pedfile    = dataset("pedigree.txt", dataset_name="demo_7animals")
genofile   = dataset("genotypes.txt", dataset_name="demo_7animals")

phenotypes = CSV.read(phenofile,DataFrame,delim=',',header=true,missingstring=["NA"])
pedigree   = get_pedigree(pedfile,separator=',',header=true)
genotypes  = get_genotypes(genofile,separator=',',method="BayesC")

# Step 2: Build model
model_equation = "y1 = intercept + x1 + x2 + x2*x3 + ID + dam + genotypes
                  y2 = intercept + x1 + x2 + ID + genotypes
                  y3 = intercept + x1 + ID + genotypes"
model = build_model(model_equation)
set_covariate(model,"x1")
set_random(model,"x2")
set_random(model,"ID dam",pedigree)

# Step 3: Define causal structure (trait1 -> trait2 and trait1 -> trait3)
my_structure = [0.0 0.0 0.0
                1.0 0.0 0.0
                1.0 0.0 0.0]

# Step 4: Run SEM-enabled MCMC
out = runMCMC(model,phenotypes; causal_structure=my_structure)
```

### Expected SEM Outputs

After `runMCMC` finishes with `causal_structure`, JWAS writes SEM-specific files:

- `structure_coefficient_MCMC_samples.txt`
- `MCMC_samples_indirect_marker_effects_genotypes_<trait>.txt`
- `MCMC_samples_overall_marker_effects_genotypes_<trait>.txt`
- `direct_marker_effects_genotypes.txt`
- `indirect_marker_effects_genotypes.txt`
- `overall_marker_effects_genotypes.txt`

`<trait>` corresponds to each phenotype in your multi-trait model (for example, `y1`, `y2`, `y3`).

### Quick Sanity Checks

- Confirm `structure_coefficient_MCMC_samples.txt` has one row per saved MCMC sample.
- Confirm indirect and overall marker-effect files exist for each trait.
- If outputs are missing, first verify the genotype term name in model equations matches generated file names.

## Reference

### `runMCMC(...; causal_structure=...)` Behavior

| Item | Behavior |
| --- | --- |
| Keyword | `causal_structure` |
| Default | `false` (SEM disabled) |
| Expected shape | square matrix with dimension = number of traits |
| Supported context | multi-trait models only |
| Value semantics | non-zero entries define directed edges (column -> row) |

### Validation and Runtime Side Effects

| Condition | Runtime behavior |
| --- | --- |
| `causal_structure` is provided in single-trait analysis | error: SEM is multi-trait only |
| `causal_structure` is not lower triangular | error: causal structure must be lower triangular |
| SEM enabled | `missing_phenotypes` is forced to `false` |
| SEM enabled | residual covariance constraint is enforced (`mme.R.constraint = true`) |

### `causal_structure` Semantics

For:

```julia
my_structure = [0.0 0.0 0.0
                1.0 0.0 0.0
                1.0 0.0 0.0]
```

- `my_structure[2,1] = 1.0` means trait 1 affects trait 2.
- `my_structure[3,1] = 1.0` means trait 1 affects trait 3.
- All entries above diagonal must be zero.

### Output Files (When SEM is Enabled)

| File | Meaning | When generated |
| --- | --- | --- |
| `structure_coefficient_MCMC_samples.txt` | sampled structural coefficient matrix entries per saved iteration | during MCMC sampling |
| `MCMC_samples_indirect_marker_effects_genotypes_<trait>.txt` | sampled indirect marker effects for trait | post-MCMC processing |
| `MCMC_samples_overall_marker_effects_genotypes_<trait>.txt` | sampled direct + indirect marker effects for trait | post-MCMC processing |
| `direct_marker_effects_genotypes.txt` | posterior summary for direct marker effects | post-MCMC processing |
| `indirect_marker_effects_genotypes.txt` | posterior summary for indirect marker effects | post-MCMC processing |
| `overall_marker_effects_genotypes.txt` | posterior summary for overall marker effects | post-MCMC processing |

### Compatibility and Limitations

| Topic | Current behavior |
| --- | --- |
| Trait count | SEM requires multi-trait setup |
| Missing phenotypes | disabled in SEM path |
| Residual covariance | constrained in SEM path |
| Marker methods | SEM post-processing expects marker-effect sample outputs for genotype terms |
| Naming dependency | output filenames include the genotype term name from model equations |
