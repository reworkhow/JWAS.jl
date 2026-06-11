# Annotated BayesC Coordinate Update Design

## Goal

Change the single-trait annotated BayesC probit coefficient update to use the same coordinate-Gibbs style as the newer nested annotation samplers: an unpenalized intercept and shrunken annotation slopes.

## Problem

After fixing the initial probit intercept, the short-chain test-data run no longer has the inflated baseline PIP. However, full annotated BayesC can still mix slowly while learning annotation effects. The current single-trait BayesC update uses one joint latent-regression Gibbs call:

```julia
rhs = ann.design_matrix' * ann.liability
Gibbs(ann.lhs, ann.coefficients, rhs, ann.variance)
```

This treats all annotation coefficients the same way and does not explicitly separate the intercept from the slopes. The nested annotation helper already uses the intended structure:

- intercept: flat prior
- annotation slopes: normal prior with variance `ann.variance`

## Proposed Fix

Update `gibbs_update_one_probit_annotation_coefficients!` to:

1. compute the latent residual `ann.liability - ann.mu`
2. call `gibbs_update_binary_probit_annotation_coefficients!`
3. rebuild `ann.mu`

The first pass will keep `ann.variance` fixed, preserving current user-facing state. Sampling an annotation-slope variance for single-trait BayesC is a larger model change and should remain separate.

## Tests

Add a unit test that sets a small annotation design with an intercept and one slope, samples liabilities from a fixed seed, and verifies the one-probit update matches `gibbs_update_binary_probit_annotation_coefficients!`. This test should fail against the old joint update and pass once single-trait BayesC uses the coordinate helper.

## Validation

Run:

```bash
julia --project=. --startup-file=no test/unit/test_annotated_bayesc.jl
```

Then rerun the short-chain `test_data` full annotated and intercept-only controls to compare with:

- before any fix
- after initialization fix only
- after initialization plus coordinate update
