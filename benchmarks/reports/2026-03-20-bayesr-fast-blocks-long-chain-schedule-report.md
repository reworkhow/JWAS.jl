# BayesR Fast Blocks Long-Chain Schedule Report

## Goal

Measure whether the production BayesR `fast_blocks` path remains acceptably close
to dense BayesR under a longer protocol, and whether a smaller explicit block size
improves that behavior.

This report compares:

- dense BayesR
- production BayesR `fast_blocks` with the current burnin-gated schedule

under:

- `chain_length = 10000`
- `burnin = 2000`
- seeds `2026:2030`
- free `pi`
- free `sigmaSq`

## Benchmark Harness

Script:

- `benchmarks/bayesr_fast_blocks_parity.jl`

Mode:

- `JWAS_BAYESR_BLOCK_MODE=long_chain_schedule_comparison`

Shared dataset:

- synthetic single-trait dense genotype matrix
- `n_obs = 60`
- `n_markers = 40`
- dataset seed `2026`
- `start_pi = [0.95, 0.03, 0.015, 0.005]`
- `gamma = [0.0, 0.01, 0.1, 1.0]`

The benchmark was run twice:

1. default `fast_blocks=true`
   - resolved block size `7`
2. explicit `fast_blocks=2`
   - resolved block size `2`

## Results

### Case 1: Default block size (`7`)

Artifacts:

- `/tmp/bayesr_fast_blocks_long_chain_default_20260320/schedule_runs.csv`
- `/tmp/bayesr_fast_blocks_long_chain_default_20260320/schedule_pairwise_summary.csv`

Burnin-gated vs dense summary:

- `sigmaSq`
  - mean abs diff: `0.043457`
  - max abs diff: `0.091404`
- residual variance
  - mean abs diff: `0.004309`
  - max abs diff: `0.009196`
- mean nonzero frequency
  - mean abs diff: `0.036223`
  - max abs diff: `0.063972`
- `pi_class1`
  - mean abs diff: `0.033362`
  - max abs diff: `0.059045`
- `pi_class2`
  - mean abs diff: `0.031128`
  - max abs diff: `0.048301`

Mean runtime:

- dense: `2.611 s`
- burnin-gated: `0.171 s`
- speedup vs dense: `15.29x`

### Case 2: Explicit block size (`2`)

Artifacts:

- `/tmp/bayesr_fast_blocks_long_chain_block2_20260320/schedule_runs.csv`
- `/tmp/bayesr_fast_blocks_long_chain_block2_20260320/schedule_pairwise_summary.csv`

Burnin-gated vs dense summary:

- `sigmaSq`
  - mean abs diff: `0.028267`
  - max abs diff: `0.057227`
- residual variance
  - mean abs diff: `0.005784`
  - max abs diff: `0.010611`
- mean nonzero frequency
  - mean abs diff: `0.022867`
  - max abs diff: `0.055212`
- `pi_class1`
  - mean abs diff: `0.021246`
  - max abs diff: `0.050424`
- `pi_class2`
  - mean abs diff: `0.022135`
  - max abs diff: `0.055263`

Mean runtime:

- dense: `2.614 s`
- burnin-gated: `0.285 s`
- speedup vs dense: `9.16x`

## Interpretation

The long-chain evidence supports three conclusions.

1. The burnin-gated BayesR block schedule is much closer to dense BayesR than
   the earlier short-chain diagnostics suggested.

2. A smaller explicit block size improves most of the BayesR summary metrics:

   - better `sigmaSq`
   - better mean nonzero frequency
   - better `pi_class1` and `pi_class2`

3. The improvement is not free:

   - `block_size = 2` is slower than the default block size
   - residual variance is slightly worse than the default block size in this run

So the current tradeoff is:

- default block size:
  - faster
  - somewhat larger occupancy and `pi` drift
- block size `2`:
  - slower
  - better agreement with dense BayesR on most BayesR-specific summaries

## Caveat About The Benchmark-Local `single_rep` Path

The same benchmark mode also recorded a benchmark-local `single_rep` path. That
path did not remain close to production dense BayesR over long chains and is not
a reliable production oracle in this form.

The useful production result in this report is therefore:

- `burnin_gated` vs `dense`

not:

- `single_rep` vs `dense`

## Conclusion

For long chains, the production BayesR burnin-gated block schedule is not showing
the catastrophic drift seen in the early short-chain diagnostics.

The remaining gap is now best understood as a speed/approximation tradeoff:

- the default block size is fastest
- a smaller block size improves agreement with dense BayesR

The next design choice is therefore policy, not basic correctness:

- keep the default block size for maximum speed
- or use a smaller BayesR block size if closer dense parity is more important
