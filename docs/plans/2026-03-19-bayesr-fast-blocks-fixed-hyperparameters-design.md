# BayesR Fast Blocks Fixed-Hyperparameters Diagnostic Design

## Goal

Diagnose the current BayesR `fast_blocks` mismatch by removing the `pi` and
`sigmaSq` feedback loops and comparing:

- dense BayesR vs BayesR `fast_blocks` with fixed `pi` and fixed `sigmaSq`
- dense BayesR vs dense BayesR across different random seeds under the same fixed
  `pi` and fixed `sigmaSq` configuration

The target is not a production feature change. The target is a clearer benchmark
that separates ordinary Monte Carlo variation from block-path-specific drift.

## Scope

This diagnostic is limited to:

- single trait
- dense genotype storage
- existing BayesR implementation
- current synthetic parity dataset

It does not add:

- a new production API
- new BayesR modeling options
- R-reference comparison for the block path

## Approach

The existing benchmark script `benchmarks/bayesr_fast_blocks_parity.jl` will be
extended with two new diagnostic modes.

### Mode 1: Fixed Hyperparameters

Purpose:

- compare dense BayesR against BayesR `fast_blocks`
- hold `pi` fixed at the current harness `start_pi`
- hold `sigmaSq` fixed at the current harness `start_sigma_sq`
- continue to sample residual variance normally

This removes two major feedback loops and makes the comparison focus on the
marker-update path and any resulting residual-variance differences.

### Mode 2: Dense Multiseed Baseline

Purpose:

- run dense BayesR only
- use the same dataset, fixed `pi`, and fixed `sigmaSq`
- vary only the random seed

This gives a direct baseline for ordinary seed-to-seed variation in BayesR.
That baseline is necessary before interpreting dense-vs-block differences as a
true block-path problem.

## Parameterization

For the fixed-hyperparameter diagnostic:

- `estimatePi = false`
- `estimate_variance = false` for BayesR marker variance
- `Pi = start_pi`
- `G_is_marker_variance = true`
- marker variance input = `start_sigma_sq`

Residual variance remains estimated in the ordinary model path.

The fixed `sigmaSq` value should be the current harness starting value, not a new
arbitrary constant, so the diagnostic stays tied to the same dataset and prior
scale already used in the parity workflow.

## Outputs

### Fixed-Hyperparameter Dense vs Block

Expected output files:

- `comparison_scalar_metrics.csv`
- `comparison_pi.csv`
- `comparison_marker_effects_top.csv`
- `runtime.csv`

Primary summaries:

- residual variance
- mean nonzero frequency
- marker-effect posterior mean differences
- marker-effect correlation

Since `pi` and `sigmaSq` are fixed, `pi.csv` should simply reflect the fixed
starting vector and `sigmaSq` should no longer be a comparison target.

### Dense Multiseed Baseline

Expected output files:

- `multiseed_runs.csv`
- `multiseed_summary.csv`

Per-seed summaries should include at least:

- residual variance
- mean nonzero frequency
- marker-effect summaries

## Minimal-Change Strategy

The diagnostic should stay in the benchmark layer.

Preferred implementation:

- extend `benchmarks/bayesr_fast_blocks_parity.jl`
- reuse `benchmarks/bayesr_parity_common.jl`
- avoid production benchmark-only flags in JWAS core code

The point is to use existing JWAS controls rather than introducing a new public
or internal sampler API solely for this investigation.

## Decision Rule

Interpretation of results:

- if dense-vs-block differences under fixed `pi` and fixed `sigmaSq` are still
  much larger than dense-vs-dense seed variation, the current block path itself
  is the likely source of mismatch
- if dense-vs-block differences shrink into the same range as dense-vs-dense
  seed variation, then the major earlier mismatch likely came from the
  `pi`/`sigmaSq` feedback loop rather than the core block update

## Deliverables

- updated `benchmarks/bayesr_fast_blocks_parity.jl`
- benchmark report:
  - `benchmarks/reports/2026-03-19-bayesr-fast-blocks-fixed-hyperparameters.md`
- implementation record:
  - `docs/plans/2026-03-19-bayesr-fast-blocks-fixed-hyperparameters-implementation.md`

## Acceptance

This diagnostic is successful if it answers the following question clearly:

Is the remaining BayesR `fast_blocks` mismatch mainly a block-kernel problem, or
is it mostly the result of the `pi` and `sigmaSq` feedback loop under block
updates?
