# BayesR Fast Blocks Fixed-Hyperparameters Diagnostic Implementation Record

## Goal

Add a diagnostic benchmark that holds BayesR `pi` and `sigmaSq` fixed so dense
vs `fast_blocks` can be compared without the mixture-weight and shared-variance
feedback loop.

## Status

Implemented and benchmarked.

This slice is diagnostic only. It does not change the production BayesR sampler.

## Files Changed

- `benchmarks/bayesr_fast_blocks_parity.jl`
- `benchmarks/bayesr_parity_common.jl`
- `test/unit/test_bayesr_parity.jl`
- `benchmarks/reports/2026-03-19-bayesr-fast-blocks-fixed-hyperparameters.md`

## What Was Added

### Benchmark Helper Support

Added helper support in `benchmarks/bayesr_parity_common.jl` for:

- fixed-hyperparameter summary writing
- runtime metadata reporting
- within-method multiseed summary aggregation

Added helper coverage in `test/unit/test_bayesr_parity.jl` for:

- fixed-hyperparameter summaries
- runtime metadata output
- within-method multiseed aggregation

### Benchmark Driver Modes

Extended `benchmarks/bayesr_fast_blocks_parity.jl` with two modes:

1. `fixed_hyperparameters`
   - dense vs block BayesR
   - fixed `pi`
   - fixed `sigmaSq`
   - residual variance still sampled

2. `dense_multiseed`
   - dense BayesR only
   - same dataset
   - fixed `pi`
   - fixed `sigmaSq`
   - multiple MCMC seeds

The driver now also:

- supports an explicit block setting through `JWAS_BAYESR_BLOCK_SETTING`
- keeps the benchmark-only block-aware burnin logic for `fast_blocks`

## Key Diagnostic Result

The earlier large block-vs-dense mismatch is not a single phenomenon.

After fixing `pi` and `sigmaSq`:

- mean nonzero frequency difference between dense and block becomes very small
  and is well within dense-only seed-to-seed variation
- residual variance difference remains materially larger than dense-only seed
  variation

This narrows the current BayesR `fast_blocks` problem considerably.

## Interpretation

The earlier evidence suggested a broad BayesR block-path mismatch.

This diagnostic shows a more precise picture:

- most of the earlier occupancy / sparsity drift was driven by the `pi` /
  `sigmaSq` feedback loop under block updates
- the remaining problem is narrower and appears downstream of that loop
- residual-variance behavior is now the main remaining mismatch signal

## Verification

Helper test command:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected benchmark commands:

```bash
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl /tmp/bayesr_fast_blocks_fixed_diag_run1

JWAS_BAYESR_BLOCK_MODE=dense_multiseed \
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl \
  /tmp/bayesr_fast_blocks_dense_multiseed_run1
```

## Recommended Follow-Up

Next diagnostic should keep `pi` and `sigmaSq` fixed and compare dense vs block:

- residual variance trace
- `yCorr` norm
- marker-effect norm

That should isolate whether the remaining discrepancy enters through:

- block exit correction
- residual-variance sampling input
- or a smaller residual difference in the marker-update path itself
