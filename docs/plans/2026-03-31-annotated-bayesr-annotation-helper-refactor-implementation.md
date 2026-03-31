# Annotated BayesR Annotation Helper Refactor Implementation

Date: 2026-03-31

## Goal

Document the two follow-up cleanup commits made after the core annotated BayesR
implementation:

- `8cccf669` `refactor: split annotation MCMC helpers`
- `04f20ae4` `refactor: clarify annotation update helpers`

These changes did not alter the supported annotated BayesR feature scope. They
were structural cleanups to make the MCMC annotation path easier to review and
maintain.

## Changes

### 1. Split annotation helper code out of `MCMC_BayesianAlphabet.jl`

Moved the annotation-specific helper block from:

- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

into a dedicated helper file:

- `src/1.JWAS/src/MCMC/annotation_updates.jl`

and added the include in:

- `src/1.JWAS/src/JWAS.jl`

Result:

- `MCMC_BayesianAlphabet.jl` now focuses on the Gibbs driver
- the annotation update logic is isolated in one file
- no sampler behavior changed in this split

### 2. Clarified and partially merged BayesC/BayesR annotation helpers

Refactored:

- `src/1.JWAS/src/MCMC/annotation_updates.jl`

Added shared binary probit primitives used by both annotated BayesC and
annotated BayesR:

- `annotation_binary_bounds!`
- `sample_binary_annotation_liabilities!`

Added clearer update boundaries:

- `gibbs_update_annotation_coefficients!(ann)` for the existing BayesC one-step
  path
- `gibbs_update_annotation_coefficients!(coeffs, X, latent_residual, variance)`
  for the BayesR stepwise scalar Gibbs updates
- `sample_annotation_effect_variance!`

Kept method-specific wrappers separate where the model semantics differ:

- BayesC one-step threshold update and prior refresh
- BayesR step-up indicator construction
- BayesR four-class `snp_pi` reconstruction

This was an intentional midpoint refactor:

- shared binary probit mechanics were merged
- method-specific logic was not over-generalized

### 3. Added equation-level documentation

Added short docstrings describing:

- binary probit latent variable interpretation
- scalar Gibbs regression updates for annotation coefficients
- BayesR step-up indicators:
  - `z1_j = 1(delta_j > 1)`
  - `z2_j = 1(delta_j > 2)`
  - `z3_j = 1(delta_j > 3)`
- BayesR reconstruction of joint class probabilities:
  - `pi_j1 = 1 - p1_j`
  - `pi_j2 = p1_j * (1 - p2_j)`
  - `pi_j3 = p1_j * p2_j * (1 - p3_j)`
  - `pi_j4 = p1_j * p2_j * p3_j`

The goal was to make the file readable as a sampler description rather than
only as code.

### 4. Small behavior-safe helper cleanup

The shared binary liability helper now accepts real-valued binary responses, not
only integer vectors. This was needed because one existing annotated BayesC test
uses `Float64` indicators.

No intended sampler behavior changed from this adjustment.

### 5. Small internal test coverage addition

Updated:

- `test/unit/test_annotated_bayesr.jl`

Added one focused internal test for:

- `annotation_binary_bounds!`

This verifies the shared truncation-bound helper independently of a full MCMC
run.

## Non-Changes

These refactors intentionally did **not** change:

- the current annotated BayesR output semantics for `pi_<genotype>`
- the current step-2 / step-3 regularization behavior
- the current BayesC latent-variance convention in the annotation sampler

Those remain separate review items.

## Verification

Executed during the refactor work:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'`
- `julia --project=. --startup-file=no test/runtests.jl`

Final full-suite result after the helper cleanup:

- `617/617` tests passed in `2m43.2s`

Docs verification for this note:

- `julia --project=docs --startup-file=no docs/make.jl`
