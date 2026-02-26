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

## Deep Technical Notes

### Core Objects and Notation

In the SEM implementation:

- `Y` is the sparse design built from observed trait values under the selected causal pattern.
- `Λ` (Lambda) is the structural coefficient matrix used to map direct and recursive relationships.
- `Λy` is the transformed phenotype vector under current SEM coefficients.

Key implementation entry points:

- `src/1.JWAS/src/structure_equation_model/SEM.jl`
- `SEM_setup`, `get_sparse_Y_FRM`, `get_Λ`, `tranform_lambda`

### Where SEM Is Invoked in MCMC

At runtime:

1. `runMCMC` stores `causal_structure` and applies SEM constraints.
2. `MCMC_BayesianAlphabet` calls `SEM_setup(...)` before sampling.
3. Each iteration, `get_Λ(...)` samples structural coefficients and updates `Λy`/corrected residual vector.
4. Saved SEM coefficient samples are written to `structure_coefficient_MCMC_samples.txt`.

Key call sites:

- `src/1.JWAS/src/JWAS.jl`
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

### Post-MCMC Marker-Effect Pipeline

After MCMC completion (SEM enabled):

1. `generate_indirect_marker_effect_sample(...)` builds per-sample indirect marker effects.
2. `generate_overall_marker_effect_sample(...)` combines direct and indirect effects.
3. `generate_marker_effect(..., effect_type)` writes posterior summaries for direct/indirect/overall effects.

These steps currently rely on consistent genotype-term naming in output marker-effect sample files.

### Current Assumptions and Caveats

- The SEM path assumes no missing-phenotype imputation during SEM sampling.
- Lower-triangular causal structure is required for identifiability/order.
- File-based post-processing is name-sensitive and depends on generated marker-effect filenames.

## Troubleshooting

### Lower-triangular matrix error

Runtime message:

`The causal structue needs to be a lower triangular matrix.`

Cause:
- `causal_structure` includes non-zero entries above the diagonal.

Fix:
- enforce lower-triangular form; keep only allowed directed edges where column affects row.

### SEM requested in single-trait model

Cause:
- SEM is enabled with only one trait.

Fix:
- use SEM only in multi-trait runs.

### Missing indirect/overall SEM output files

Cause:
- genotype term name in model equations does not match expected marker-effect file naming.

Fix:
- keep genotype term naming consistent (for example, `genotypes` in both model equation and output expectations).

### Historical regression (issue #162)

`runMCMC(...; causal_structure=...)` previously failed from an internal field mismatch.
This path is covered by unit regression test `test/unit/test_sem_issue162.jl`.

## Related Pages

- [Workflow](workflow.md)
- [Public API](public.md)
- [Block BayesC](block_bayesc.md)
