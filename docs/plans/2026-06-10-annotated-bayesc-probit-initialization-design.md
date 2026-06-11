# Annotated BayesC Probit Initialization Design

## Goal

Make single-trait annotated BayesC start its annotation probit state from the user-supplied BayesC `Pi` instead of the implicit `Phi(0) = 0.5` inclusion probability.

## Problem

Single-trait annotated BayesC currently preserves the starting marker exclusion probabilities in `genotypei.π`, but `MarkerAnnotations` initializes annotation coefficients and `mu` to zero. For sparse BayesC settings such as `Pi = 0.99`, this creates an inconsistent starting state:

- marker sampler starts from inclusion probability `0.01`
- annotation probit model starts from inclusion probability `0.5`

The test-data diagnostics showed this produces slow burn-in in intercept-only annotated BayesC. With a short chain, the baseline PIP averaged around `0.05`; with a longer chain, it converged back near the plain BayesC level around `0.01`.

## Proposed Fix

In `initialize_bayesc_single_trait_annotations!`, after rebuilding `MarkerAnnotations`, initialize the annotation intercept to match the starting inclusion probability:

```julia
start_inclusion = mean(1 .- start_pi)
intercept = quantile(Normal(), clamp(start_inclusion, eps(Float64), 1 - eps(Float64)))
ann.coefficients[1] = intercept
ann.mu .= ann.design_matrix * ann.coefficients
```

Annotation slopes remain zero. For marker-specific starting `Pi` vectors, a single intercept cannot match every marker, so the intercept matches the average starting inclusion probability.

## Tests

Update the annotated BayesC startup unit test to verify:

- `genotypei.π` still stores the requested starting exclusion probabilities
- the annotation intercept is `qnorm(mean(1 - Pi))`
- annotation slopes remain zero
- `ann.mu` is initialized from `design_matrix * coefficients`
- repeated initialization is deterministic

## Validation

Run the focused unit test file:

```bash
julia --project=. --startup-file=no test/unit/test_annotated_bayesc.jl
```

Then rerun the single-trait test-data benchmark to check the impact on short-chain annotated BayesC behavior.
