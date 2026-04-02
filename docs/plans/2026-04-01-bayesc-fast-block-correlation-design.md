# BayesC Fast-Block Correlation Benchmark Design

## Goal

Measure how closely `fast_blocks=1` tracks the dense production sampler for:

- standard `BayesC`
- annotated `BayesC`

using the same synthetic dataset and the same MCMC seed for each dense/block
pair.

## Question

The user expectation is:

- dense `BayesC` and `BayesC` with `fast_blocks=1` should be extremely close
- dense annotated `BayesC` and annotated `BayesC` with `fast_blocks=1` should
  also remain highly correlated even if they are not bitwise identical

The benchmark should quantify that closeness directly instead of relying on
qualitative impressions.

## Approach

Use the production `runMCMC` path on a moderate synthetic single-trait dataset.

For each seed:

1. run dense `BayesC`
2. run `BayesC` with `fast_blocks=1`
3. run dense annotated `BayesC`
4. run annotated `BayesC` with `fast_blocks=1`

All four runs use:

- the same genotype matrix
- the same phenotype vector
- the same chain settings
- paired same-seed comparisons within method

## Metrics

For each dense/block pair, compute:

- correlation of posterior mean marker-effect estimates
- correlation of posterior model frequencies
- correlation of EBV predictions
- residual-variance absolute difference
- BayesC `pi` absolute difference for standard BayesC
- per-marker annotation-induced `pi_j` correlation for annotated BayesC
- annotation-coefficient correlation and max absolute difference for annotated
  BayesC

## Dataset

Use one synthetic annotated benchmark dataset with:

- `n_obs = 200`
- `n_markers = 1000`
- `h2 = 0.45`
- scenario `stepwise_annotation_signal`

This gives a realistic annotation signal while still being small enough to run
paired dense/block comparisons quickly.

## Deliverables

- benchmark script:
  - `benchmarks/bayesc_fast_block_correlation.jl`
- report:
  - `benchmarks/reports/2026-04-01-bayesc-fast-block-correlation-report.md`
- implementation note:
  - `docs/plans/2026-04-01-bayesc-fast-block-correlation-implementation.md`
