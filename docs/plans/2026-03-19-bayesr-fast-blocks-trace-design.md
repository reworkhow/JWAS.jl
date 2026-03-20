# BayesR Fast Blocks Trace Design

## Goal

Add a short diagnostic trace for BayesR `fast_blocks` under:

- fixed `pi`
- fixed `sigmaSq`
- estimated residual variance

The purpose is to isolate where the remaining dense-vs-block mismatch first
appears after the earlier fixed-hyperparameter diagnostic narrowed it to
residual-variance-related behavior.

## Scope

This is benchmark-only.

It does not change the production JWAS BayesR sampler or public API.

The trace benchmark will:

- use the existing synthetic single-trait dense dataset path
- compare dense BayesR vs BayesR `fast_blocks`
- use the same MCMC seed in both runs
- run a short chain:
  - `chain_length = 100`
  - `burnin = 0`

## Benchmark Mode

Extend `benchmarks/bayesr_fast_blocks_parity.jl` with a new mode:

- `JWAS_BAYESR_BLOCK_MODE=fixed_hyperparameters_trace`

This mode will:

- keep `pi` fixed
- keep `sigmaSq` fixed
- still sample residual variance
- write per-iteration trace files for dense and block BayesR
- write one comparison file aligned by iteration

## Trace Contents

Per iteration, record:

- `iter`
- `residual_variance`
- `ycorr_norm`
- `alpha_norm`
- `alpha_abs_mean`
- `nnz`
- `max_abs_alpha`

These are the minimum diagnostics needed to decide whether the remaining block
drift enters through:

- upstream state propagation (`yCorr`, `alpha`)
- or downstream residual-variance behavior

## Implementation Approach

Keep the trace harness benchmark-local.

Do not modify the production `runMCMC` path to expose per-iteration internals.

Instead:

- reuse the existing synthetic dataset and fixed-hyperparameter setup
- add a lightweight local Gibbs loop in the benchmark script
- call the production BayesR marker kernels:
  - `BayesR!`
  - `BayesR_block!`
- after each iteration:
  - update residual variance locally
  - record the trace row

This keeps the implementation small while still exercising the production dense
and block marker-update kernels.

## Decision Rule

Interpret the trace in this order:

1. if `ycorr_norm` diverges before `residual_variance`, the issue is upstream in
   state propagation
2. if `alpha` summaries stay aligned but `residual_variance` diverges, the issue
   is concentrated in the residual-variance update inputs
3. if `nnz` stays aligned, the remaining mismatch is not an occupancy problem in
   this fixed-hyperparameter setting

## Deliverables

- benchmark update:
  - `benchmarks/bayesr_fast_blocks_parity.jl`
- trace report:
  - `benchmarks/reports/2026-03-19-bayesr-fast-blocks-trace.md`
- implementation record:
  - `docs/plans/2026-03-19-bayesr-fast-blocks-trace-implementation.md`
