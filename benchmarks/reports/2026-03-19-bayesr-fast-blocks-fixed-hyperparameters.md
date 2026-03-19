# BayesR Fast Blocks Fixed-Hyperparameters Diagnostic

## Goal

Re-run the BayesR dense-vs-`fast_blocks` comparison with:

- fixed `pi`
- fixed `sigmaSq`
- estimated residual variance

Then compare that block-vs-dense difference against ordinary dense-BayesR
seed-to-seed variation on the same dataset.

The purpose is to separate:

- block-kernel drift
from
- drift caused by the `pi` / `sigmaSq` feedback loop

## Setup

Script:

- `benchmarks/bayesr_fast_blocks_parity.jl`

Mode 1:

- `JWAS_BAYESR_BLOCK_MODE=fixed_hyperparameters`

Mode 2:

- `JWAS_BAYESR_BLOCK_MODE=dense_multiseed`

Shared dataset:

- synthetic single-trait dense dataset
- `n_obs = 60`
- `n_markers = 40`
- dataset seed `2026`
- `start_pi = [0.95, 0.03, 0.015, 0.005]`
- `gamma = [0.0, 0.01, 0.1, 1.0]`
- `start_h2 = 0.5`

Fixed-hyperparameter settings:

- `estimatePi = false`
- BayesR marker variance estimation disabled
- marker variance fixed at `start_sigma_sq`
- residual variance still sampled

## Results

### Dense vs Fast Blocks With Fixed `pi` And Fixed `sigmaSq`

Command:

```bash
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl /tmp/bayesr_fast_blocks_fixed_diag_run1
```

Block configuration:

- `fast_blocks=true`
- resolved block size: `7`
- requested burnin: `100`
- effective block burnin in benchmark harness: `14`

Runtime:

- dense: `9.0412 s`
- fast_blocks: `0.4778 s`
- speedup: `18.92x`

Posterior-summary comparison:

- residual variance
  - dense: `0.960822`
  - block: `0.922044`
  - abs diff: `0.038778`
  - rel diff: `4.04%`
- mean nonzero frequency
  - dense: `0.047875`
  - block: `0.048214`
  - abs diff: `0.000339`
  - rel diff: `0.71%`
- fixed `pi`
  - identical in both paths

Interpretation:

- the large earlier nonzero-frequency drift is gone once `pi` and `sigmaSq` are
  fixed
- the remaining visible mismatch is mostly in residual variance

### Dense BayesR Seed-to-Seed Variation

Command:

```bash
JWAS_BAYESR_BLOCK_MODE=dense_multiseed \
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl \
  /tmp/bayesr_fast_blocks_dense_multiseed_run1
```

Seeds:

- `2026, 2027, 2028, 2029, 2030`

Observed dense-only runs:

- residual variance range:
  - min: `0.950336`
  - max: `0.964192`
  - range: `0.013856`
- mean nonzero frequency range:
  - min: `0.041375`
  - max: `0.047875`
  - range: `0.006500`
- mean absolute marker effect range:
  - min: `0.004331`
  - max: `0.005136`
  - range: `0.000806`

Note:

- the `seconds` column in the multiseed output is not useful for inference
  comparison because the first run includes JIT compilation cost

## Comparison Against Seed Variation

Fixed-hyperparameter dense-vs-block differences vs dense-only seed spread:

- residual variance
  - block-vs-dense abs diff: `0.038778`
  - dense seed range: `0.013856`
  - block difference is materially larger than ordinary seed variation
- mean nonzero frequency
  - block-vs-dense abs diff: `0.000339`
  - dense seed range: `0.006500`
  - block difference is well inside ordinary seed variation

## Conclusion

This diagnostic changes the interpretation of the earlier BayesR `fast_blocks`
mismatch.

What it shows:

- the large earlier drift in nonzero frequency and class occupancy was mostly
  caused by the `pi` / `sigmaSq` feedback loop under block updates
- once `pi` and `sigmaSq` are fixed, the block path closely matches dense BayesR
  on marker occupancy
- a residual-variance mismatch still remains, and it is larger than ordinary
  seed-to-seed variation in dense BayesR

So the current BayesR `fast_blocks` issue is now much narrower:

- not a broad failure of the block marker update
- likely a remaining interaction between the block path and residual-variance
  behavior, or another downstream state difference after the marker sweep

## Recommended Next Step

Keep `pi` and `sigmaSq` fixed and trace:

- residual variance
- `yCorr` norm
- marker-effect norm

for dense vs block BayesR on the same dataset and seed.

That should isolate whether the remaining mismatch comes from:

- block exit correction to `yCorr`
- residual-variance sampling input
- or a smaller remaining discrepancy in block marker updates
