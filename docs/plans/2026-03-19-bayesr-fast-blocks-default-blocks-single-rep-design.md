# BayesR Fast Blocks Default-Blocks Single-Rep Design

## Goal

Add a benchmark-only diagnostic that keeps the default BayesR block partitioning but forces one within-block sweep per outer iteration. The diagnostic should compare this block-cached path against dense BayesR using the same outer chain length, free `pi`, and free `sigmaSq`.

## Motivation

Current BayesR fast-block behavior still differs from dense BayesR for `block_size > 1`. We already established that:

- `fast_blocks = 1` matches dense BayesR essentially exactly
- the remaining drift enters when true multi-marker block behavior is used

This diagnostic isolates whether the problem comes from repeated within-block sweeps (`nreps = block_size`) or from the default block partition / block-exit correction itself.

## Scope

- Benchmark-only change
- Default block size (`floor(sqrt(nObs))`)
- One within-block sweep per outer iteration (`nreps = 1`)
- Same outer `chain_length` and `burnin` as dense BayesR
- Free `pi`
- Free `sigmaSq`
- Free residual variance
- Single-trait only

No production sampler behavior changes are included in this slice.

## Approach Options

### Option 1: Benchmark-local Gibbs driver

Implement a local Gibbs loop in `benchmarks/bayesr_fast_blocks_parity.jl` for:

- dense BayesR
- default-blocks BayesR with one within-block sweep

Reuse production math where possible:

- `BayesR!`
- `BayesR_block!`
- `samplePi`
- `sample_marker_effect_variance`
- `sample_variance`

Recommendation: use this option. It keeps the diagnostic local and compares only the marker-sweep schedule.

### Option 2: Temporary production switch

Add a temporary production flag to force `nreps = 1`.

Rejected: too much surface area for a diagnostic.

## Design

Add a new benchmark mode to `benchmarks/bayesr_fast_blocks_parity.jl`:

- `default_blocks_single_rep`

The mode should:

1. Build the usual synthetic BayesR dataset.
2. Initialize two fresh local BayesR states:
   - dense
   - default-block partition
3. Run the same local outer Gibbs loop for both cases:
   - update `mu`
   - update marker effects
   - update `pi`
   - update `sigmaSq`
   - update residual variance
4. For the block case, call the current production `BayesR_block!` overload with a very large `burnin` so `bayesr_block_nreps(iter, burnin, block_size)` stays at `1` for all iterations.
5. Accumulate posterior summaries after burnin and write comparison CSVs.

## Outputs

The mode should write:

- `comparison_scalar_metrics.csv`
- `comparison_pi.csv`
- `comparison_marker_effects_top.csv`
- `runtime.csv`

## Decision Rule

- If default-blocks single-rep matches dense BayesR closely, then the main remaining drift is caused by repeated within-block sweeps, not by the default block partition itself.
- If it still drifts materially, then the default block partition / block-exit correction is still changing the chain even with one sweep per SNP.
