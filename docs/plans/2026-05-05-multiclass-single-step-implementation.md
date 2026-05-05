# Multi-Class Single-Step Implementation

## Summary

Implemented minimal dense multi-genotype-category support for single-step BayesC.
The old one-category path is preserved, including internal `ϵ` and `J` term names.

## Changes

- Removed the validation error that rejected multiple genotype categories in
  single-step analysis.
- Updated `SSBRrun` to impute each `Genotypes` object in `mme.M`.
- Kept one shared `ϵ` and one shared `J` for multi-category single-step runs.
- Set the shared `ϵ` genetic variance to the sum of the genotype categories'
  genomic variance values.
- Added a regression test for `y1 = intercept + geno1 + geno2` with
  `single_step_analysis=true`.

## Verification

- `julia --project=. --startup-file=no -e 'using Test; include("test/unit/test_single_step.jl")'`
- `julia --project=. --startup-file=no -e 'using Test; include("test/unit/test_input_validation.jl")'`
- `julia --project=. --startup-file=no -e 'using Test; include("test/unit/test_set_random.jl")'`
- `git diff --check`
