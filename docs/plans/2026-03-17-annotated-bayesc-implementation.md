# Annotated BayesC Implementation

**Date:** 2026-03-17

## Goal

Document how Annotated BayesC was implemented in JWAS on top of the current single-trait BayesC code path, including the final architecture, MCMC update order, storage-mode support, and the implementation decisions that were settled during review.

## Scope

### Implemented
- Single-trait `method="BayesC"` support for annotation-driven prior exclusion probabilities.
- Dense genotype support.
- Dense `fast_blocks` support.
- `storage=:stream` support.
- Annotation-specific validation, setup, posterior updates, output, and tests.
- Backward-compatibility adjustments so ordinary BayesC behavior stays stable for non-annotated runs.

### Not implemented
- Multi-trait Annotated BayesC.
- Annotated methods other than BayesC.
- `storage=:stream` together with `fast_blocks`.
- Annotation priors for GBLUP, BayesA, BayesB, BayesL, or BayesR.

## User-Facing API

Annotated BayesC is enabled through:

```julia
get_genotypes(...; method="BayesC", annotations=annotation_matrix)
```

Key interface rules:
- `annotations` is a numeric matrix with one row per raw marker.
- Users do not supply an intercept column.
- JWAS prepends the intercept internally after marker QC/filtering.
- Annotated runs require `estimatePi=true`. If users pass `estimatePi=false`, JWAS warns and overrides it.

Relevant code:
- `src/1.JWAS/src/markers/readgenotypes.jl`

## Internal Representation

The implementation keeps annotation state separate from the core genotype matrix by storing a `MarkerAnnotations` object on the `Genotypes` struct.

That internal state holds:
- filtered annotation design matrix with the prepended intercept
- current annotation coefficients
- marker liabilities for the threshold model
- lower and upper truncation bounds
- linear predictor `mu`
- liability variance
- precomputed `X'X`

This avoids spreading annotation-specific arrays across `Genotypes` and keeps Annotated BayesC state localized.

Relevant code:
- `src/1.JWAS/src/types.jl`

## Prior Representation

Standard single-trait BayesC conceptually uses one exclusion probability `pi`. Annotated BayesC needs marker-specific exclusion probabilities `pi_j`.

The implementation resolves that by storing single-trait BayesC `pi` internally as a length-`nMarkers` vector:
- ordinary BayesC uses `fill(pi_scalar, nMarkers)`
- Annotated BayesC updates the vector marker-by-marker

This matches the existing pattern used for marker-effect variances, where BayesC can use repeated values while BayesB already uses true per-marker vectors.

Important consequence:
- dense, block, and streaming BayesC paths all consume the same per-marker `pi[j]` representation
- output and screen-print logic had to be adjusted so non-annotated BayesC still looks like scalar `pi` to users

Relevant code:
- `src/1.JWAS/src/markers/tools4genotypes.jl`
- `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`
- `src/1.JWAS/src/output.jl`
- `src/1.JWAS/src/JWAS.jl`

## Statistical Form

Standard BayesC uses:
- `P(delta_j = 0) = pi`
- `P(delta_j = 1) = 1 - pi`

Annotated BayesC replaces the common exclusion probability with a probit threshold model:

- `z_j = x_j' * gamma + e_j`
- `e_j ~ N(0, sigma_ann^2)`
- `delta_j = 1` if `z_j > 0`, else `0`
- `pi_j = P(delta_j = 0 | x_j, gamma) = 1 - Phi((x_j' * gamma) / sigma_ann)`

So annotations do not enter the phenotype likelihood directly. They update the prior odds that a marker is excluded.

## Input Validation And Setup

Annotation validation happens in `get_genotypes`, before MCMC starts.

### Raw validation
- `annotations` must be false or a numeric matrix
- `method` must be `"BayesC"`
- annotation row count must match the number of raw markers

### Design validation
After QC/filtering and before building `MarkerAnnotations`, JWAS rejects:
- constant user annotation columns
- exact collinearity after the intercept is added

This was added because the annotation coefficient update uses `X'X` directly and does not impose extra shrinkage on annotation coefficients.

### Dense path
Dense loading reads the full marker matrix, performs QC, applies the same marker mask to `annotations`, prepends the intercept, and builds `MarkerAnnotations`.

### Streaming path
Streaming support required backend metadata to preserve:
- total raw marker count before QC
- retained raw-marker indices after QC

That mapping is needed so raw-marker annotations can be filtered exactly like dense mode. Annotated streaming now errors clearly if a legacy backend without this metadata is used.

Relevant code:
- `src/1.JWAS/src/markers/readgenotypes.jl`
- `src/1.JWAS/src/markers/streaming_genotypes.jl`

## MCMC Update Order

Annotated BayesC is integrated into the existing `MCMC_BayesianAlphabet` driver rather than using a separate sampler.

The update order is:

1. Initialize marker inclusion indicators `delta` for annotated runs.
2. Perform one annotation-prior update before the first main BayesC sweep.
3. On each MCMC iteration:
   - update phenotype liabilities first if the phenotype is categorical/censored
   - update marker effects and `delta` through the normal BayesC marker step
   - update annotation liabilities, annotation coefficients, and per-marker `pi_j`
4. Use the refreshed `pi_j` vector on the next BayesC sweep.

This is an alternating Gibbs scheme where:
- phenotype data update `delta`
- `delta` updates the annotation threshold model
- the annotation model updates `pi_j`
- `pi_j` feeds back into the next marker update

Relevant code:
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

## Annotation Initialization

Initialization was one of the main design points settled during review.

### Final rule
- If starting `pi == 0.0`, initialize with 10% excluded markers and 90% included markers.
- If starting `pi == 1.0`, initialize with 10% included markers and 90% excluded markers.
- If `0 < pi < 1`, sample `delta_j` independently from the starting exclusion probability.
- If that random draw still produces all included or all excluded markers, repair it with the same 10% minority rule.

JWAS prints an `@info` message whenever one of these degenerate-start fallbacks is used.

### Why
The annotation threshold model cannot start from a degenerate `delta` vector where every marker is in the same class.

Relevant code:
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

## Annotation Liability And Coefficient Updates

The annotation block performs three steps:

1. Convert current `delta` values into truncation bounds.
2. Sample marker liabilities from truncated normals.
3. Regress liabilities on the annotation design matrix to update annotation coefficients.

Then:
- `mu = X * coefficients`
- `pi_j = 1 - Phi(mu_j / sqrt(variance))`

The final semantics are:
- `ann.variance` is a variance
- `Normal(...)` calls use `sqrt(ann.variance)`
- `Gibbs(...)` receives `ann.variance` directly

This was corrected explicitly during review to keep naming and usage consistent.

Relevant code:
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

## BayesABC Integration

The existing BayesABC marker-update code was generalized so all single-trait execution paths use per-marker `pi[j]`:
- dense
- `fast_blocks`
- streaming

The shared helper `bayesabc_update_marker!` receives:
- marker-specific exclusion log-probability `log(pi_j)`
- marker-specific inclusion log-probability `log(1 - pi_j)`

To keep the shared code safe, `bayesabc_pi_vector` normalizes scalar-vs-vector `pi` and errors clearly if a provided `pi` vector length does not match the number of markers.

Relevant code:
- `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`

## Compatibility Decisions

Several implementation details were adjusted to preserve stable behavior for existing non-annotated BayesC runs.

### `estimatePi`
- Annotated BayesC forces `estimatePi=true`
- if users pass `estimatePi=false`, JWAS warns and overrides it

### Non-annotated output
Even though ordinary single-trait BayesC now uses a vector internally, user-facing output for non-annotated runs remains scalar:
- summary output collapses back to one `pi`
- MCMC sample files for ordinary BayesC still behave like scalar-`pi` output

### Annotated screen print
Annotated runs do not print the full `pi_j` vector to screen. Instead JWAS prints:

- `pi_j (min/mean/max)`

That keeps the output readable for large marker sets while still acknowledging the per-marker state.

Relevant code:
- `src/1.JWAS/src/markers/readgenotypes.jl`
- `src/1.JWAS/src/output.jl`
- `src/1.JWAS/src/JWAS.jl`

## Streaming-Specific Decision

Annotated `storage=:stream` now requires a backend prepared with the current `prepare_streaming_genotypes`, because the backend must include raw-marker mapping metadata.

JWAS does not silently guess compatibility for old backends. Instead it errors clearly and tells the user to rebuild the backend.

This decision was made to avoid silent misalignment between streamed markers and the annotation matrix.

Relevant code:
- `src/1.JWAS/src/markers/streaming_genotypes.jl`
- `src/1.JWAS/src/markers/readgenotypes.jl`

## Output And Reporting

Annotated BayesC adds posterior summaries for annotation coefficients. The user doc focuses on the method and how to run it, while the implementation keeps developer-facing annotation bookkeeping in the output layer.

Output behavior:
- annotation coefficients are accumulated with posterior mean and second-moment tracking
- ordinary BayesC retains scalar-looking `pi`
- annotated runs keep per-marker `pi_j` in output files

Relevant code:
- `src/1.JWAS/src/output.jl`

## Testing Strategy

Annotated BayesC added focused unit coverage for:
- API validation
- unsupported combinations
- annotation input shape checks
- constant and collinear annotation design rejection
- dense runs
- `fast_blocks` runs
- streaming runs
- streaming backend compatibility
- `pi` vector length guard
- output and print behavior
- initialization and annotation-prior update details

Relevant code:
- `test/unit/test_annotated_bayesc.jl`
- `test/runtests.jl`

## Files Touched By The Implementation

Core implementation:
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`
- `src/1.JWAS/src/markers/readgenotypes.jl`
- `src/1.JWAS/src/markers/streaming_genotypes.jl`
- `src/1.JWAS/src/types.jl`
- `src/1.JWAS/src/output.jl`
- `src/1.JWAS/src/JWAS.jl`

Tests and docs:
- `test/unit/test_annotated_bayesc.jl`
- `docs/src/manual/annotated_bayesc.md`
- `docs/src/manual/workflow.md`

## Final Implementation Summary

Annotated BayesC was implemented as an extension of the current single-trait BayesC pipeline, not as a forked sampler. The key architectural choice was to represent BayesC `pi` internally as a per-marker vector, then let annotations update that vector through a marker-level probit threshold model inside the existing MCMC flow.

Dense, block, and streaming BayesC now all consume the same per-marker prior representation. Input validation, backward-compatibility behavior, and streaming metadata checks were added so that annotation support integrates cleanly without regressing ordinary BayesC behavior.
