# BayesR Fast Blocks Parity Report

## Goal

Compare single-trait dense BayesR against the new single-trait dense `fast_blocks`
BayesR path in JWAS. The acceptance target for this first cut was posterior-summary
parity against dense BayesR, not parity against the external R reference.

## Benchmark Harness

Script:

- `benchmarks/bayesr_fast_blocks_parity.jl`

Dataset:

- synthetic dense genotype matrix
- `n_obs = 60`
- `n_markers = 40`
- one trait
- starting `pi = [0.95, 0.03, 0.015, 0.005]`
- `gamma = [0.0, 0.01, 0.1, 1.0]`

Compared summaries:

- posterior mean `sigmaSq`
- posterior mean residual variance
- posterior mean nonzero marker frequency
- posterior mean `pi`
- marker-effect posterior means and nonzero frequencies
- runtime

## Important Benchmark Detail

`runMCMC(...; fast_blocks=...)` rescales `chain_length` internally by block size.
It does **not** rescale `burnin`.

For the parity harness, the block run therefore uses a block-aware burnin to avoid
the degenerate case where `burnin` exceeds the internally rescaled outer chain
length and no post-burnin samples are written.

This is a benchmarking harness adjustment only. It does not change the production
JWAS `fast_blocks` behavior.

## Results

### Case 1: Default `fast_blocks=true`, short chain

Command:

```bash
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl
```

Observed block size:

- `7`

Runtime summary:

- dense: `9.5433 s`
- fast_blocks: `1.6017 s`
- speedup: `5.96x`

Posterior-summary comparison:

- `sigmaSq` relative difference: `90.9%`
- residual variance relative difference: `3.77%`
- mean nonzero frequency relative difference: `79.2%`
- max `pi` abs difference: `0.3436`

Conclusion:

- not acceptable
- default block sizing is too aggressive in this setting

### Case 2: Default `fast_blocks=true`, longer chain

Command:

```bash
JWAS_BAYESR_BLOCK_CHAIN_LENGTH=2100 \
JWAS_BAYESR_BLOCK_BURNIN=700 \
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl \
  /tmp/bayesr_fast_blocks_parity_long
```

Observed block size:

- `7`

Runtime summary:

- dense: `9.6815 s`
- fast_blocks: `0.5942 s`
- speedup: `16.29x`

Posterior-summary comparison:

- `sigmaSq` relative difference: `40.5%`
- residual variance relative difference: `1.81%`
- mean nonzero frequency relative difference: `21.9%`
- class-1 `pi` abs difference: `0.1205`
- class-3 `pi` abs difference: `0.1299`

Conclusion:

- longer chains reduce the gap
- parity is still not acceptable

### Case 3: Smaller explicit block size, longer chain

Command:

```bash
JWAS_BAYESR_BLOCK_CHAIN_LENGTH=2100 \
JWAS_BAYESR_BLOCK_BURNIN=700 \
JWAS_BAYESR_BLOCK_SETTING=2 \
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl \
  /tmp/bayesr_fast_blocks_parity_block2
```

Observed block size:

- `2`

Runtime summary:

- dense: `9.6986 s`
- fast_blocks: `0.6293 s`
- speedup: `15.41x`

Posterior-summary comparison:

- `sigmaSq` relative difference: `29.5%`
- residual variance relative difference: `2.30%`
- mean nonzero frequency relative difference: `17.3%`
- class-1 `pi` abs difference: `0.0895`
- class-2 `pi` abs difference: `0.0674`

Conclusion:

- smaller blocks materially improve parity
- the current BayesR block path is still not close enough to dense BayesR

## Interpretation

What is clear from the current evidence:

- the benchmark harness itself now works correctly
- the earlier all-zero block summaries were a harness bug caused by unscaled
  burnin against internally rescaled `chain_length`
- after fixing the harness, the BayesR block path still shows material posterior
  drift relative to dense BayesR
- the drift is strongest in `sigmaSq`, class occupancy, and nonzero frequency
- smaller explicit block sizes help, which points to a block-schedule or
  block-kernel sensitivity rather than a reporting issue

## Current Status

The current BayesR `fast_blocks` implementation is functional and fast, but it
does **not** yet meet the planned posterior-summary parity acceptance criterion.

It should not be treated as merge-ready on the basis of the current benchmark
evidence.

## Recommended Next Debug Step

1. Run dense-vs-block BayesR with `estimatePi=false` to isolate the block marker
   update from the Dirichlet `pi` update.
2. Compare dense vs block trajectories for:
   - `sigmaSq`
   - `nnz`
   - class counts
   - `pi`
3. Revisit the BayesR block schedule rather than assuming the BayesC schedule is
   valid unchanged.
