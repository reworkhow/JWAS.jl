# BayesC and BayesR Fast-Blocks Comparison Report

## Goal

Compare the production dense and `fast_blocks` paths for single-trait BayesC and
BayesR under the same long-chain multiseed protocol.

## Protocol

- benchmark script:
  [bayesc_bayesr_fast_blocks_comparison.jl](/Users/haocheng/Github/JWAS.jl/benchmarks/bayesc_bayesr_fast_blocks_comparison.jl)
- dataset:
  synthetic single-trait data from
  [bayesr_parity_common.jl](/Users/haocheng/Github/JWAS.jl/benchmarks/bayesr_parity_common.jl)
- `n_obs = 60`
- `n_markers = 40`
- `chain_length = 10000`
- `burnin = 2000`
- seeds:
  `2026,2027,2028,2029,2030`

Cases:

- `BayesC_dense`
- `BayesC_fast_blocks_default`
- `BayesC_fast_blocks_1`
- `BayesR_dense`
- `BayesR_fast_blocks_default`
- `BayesR_fast_blocks_1`

Output directory used for the real run:

- `/tmp/bayesc_bayesr_fast_blocks_comparison_ebv_20260321`

## Results

### BayesC: `fast_blocks=1` vs dense

Mean absolute differences across seeds:

- residual variance: `4.56e-6`
- mean model frequency: `0.0`
- marker abs mean: `5.14e-8`
- phenotype-EBV correlation: `2.51e-7`
- scalar `pi`: `0.0`
- marker variance: `3.12e-8`

Interpretation:

- for BayesC, `fast_blocks=1` is effectively identical to dense BayesC

### BayesC: default `fast_blocks=true` vs dense

Mean absolute differences across seeds:

- residual variance: `9.14e-3`
- mean model frequency: `4.41e-2`
- marker abs mean: `2.56e-3`
- phenotype-EBV correlation: `1.03e-2`
- scalar `pi`: `4.34e-2`
- marker variance: `1.74e-3`

Interpretation:

- default BayesC `fast_blocks` changes the posterior summaries materially relative
  to dense BayesC

### BayesR: `fast_blocks=1` vs dense

Mean absolute differences across seeds:

- residual variance: `1.24e-3`
- mean model frequency: `4.12e-3`
- marker abs mean: `2.59e-4`
- phenotype-EBV correlation: `3.18e-4`
- `sigmaSq`: `4.24e-3`
- `pi_class1`: `4.11e-3`
- `pi_class2`: `4.09e-3`
- `pi_class3`: `7.37e-4`
- `pi_class4`: `2.15e-4`

Interpretation:

- `fast_blocks=1` is much closer to dense BayesR than the default block path
- it is not numerically exact under free hyperparameter updates, but it is a
  close approximation

### BayesR: default `fast_blocks=true` vs dense

Mean absolute differences across seeds:

- residual variance: `4.31e-3`
- mean model frequency: `3.62e-2`
- marker abs mean: `1.48e-3`
- phenotype-EBV correlation: `3.14e-3`
- `sigmaSq`: `4.35e-2`
- `pi_class1`: `3.34e-2`
- `pi_class2`: `3.11e-2`
- `pi_class3`: `9.08e-3`
- `pi_class4`: `2.96e-3`

Interpretation:

- default BayesR `fast_blocks` also changes posterior summaries materially
- the largest BayesR shifts are in `sigmaSq` and class probabilities

## Prediction Accuracy

Mean in-sample `cor(y, EBV)` across seeds:

- `BayesC_dense`: `0.665850`
- `BayesC_fast_blocks_default`: `0.665609`
- `BayesC_fast_blocks_1`: `0.665850`
- `BayesR_dense`: `0.676784`
- `BayesR_fast_blocks_default`: `0.677397`
- `BayesR_fast_blocks_1`: `0.677102`

Interpretation:

- for BayesC, `fast_blocks=1` is effectively identical to dense on phenotype-EBV correlation
- default BayesC `fast_blocks` moves the correlation slightly more, but still by a small amount in absolute terms
- for BayesR, both block variants stay very close to dense on phenotype-EBV correlation
- the prediction-correlation differences are smaller than the BayesR shifts seen in `sigmaSq` and class probabilities

## Runtime Note

The raw runtime columns were saved in `comparison_runs.csv`, but these timings
include JIT warmup effects because all cases ran in one Julia process. They are
useful for rough context only and should not be treated as authoritative speed
measurements.

## Conclusions

1. `fast_blocks=1` is effectively dense-equivalent for BayesC.
2. `fast_blocks=1` is also close to dense for BayesR, but not exact under the
   full free-hyperparameter production protocol.
3. Default `fast_blocks=true` changes both methods relative to dense.
4. BayesR remains more sensitive than BayesC in the hyperparameter summaries
   that matter most to BayesR, especially `sigmaSq` and class probabilities.
5. On in-sample phenotype-EBV correlation, both methods are more stable than the
   variance and mixture summaries, especially for BayesR.

## Files

- runs:
  `/tmp/bayesc_bayesr_fast_blocks_comparison_ebv_20260321/comparison_runs.csv`
- pairwise summary:
  `/tmp/bayesc_bayesr_fast_blocks_comparison_ebv_20260321/comparison_pairwise_summary.csv`
