# Multi-Class Single-Step Design

## Goal

Allow single-step BayesC runs to include more than one dense genotype category,
for example `y = intercept + geno1 + geno2`, without changing the existing
single-genotype-category path.

## Minimal Approach

- Keep the existing dense genotype requirements, including the current rule that
  all genotype categories must have the same individual IDs.
- Remove only the single-step validation gate that rejects `length(mme.M) != 1`.
- Reuse the existing pedigree inverse and imputation code for each genotype
  category in `mme.M`.
- Preserve the current internal terms `ϵ` and `J` for both single-category and
  multi-category analyses.
- Use one shared `ϵ` with genetic variance equal to the sum of the genotype
  categories' genomic variance values.
- Keep streaming, GBLUP single-step, and different-ID genotype categories out of
  scope.

## Verification

Add a unit regression test that builds two dense BayesC genotype categories from
the demo data and runs `single_step_analysis=true`. The test must fail on the
current validation error before implementation and pass after the minimal change.
