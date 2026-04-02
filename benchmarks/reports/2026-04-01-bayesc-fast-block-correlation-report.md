# BayesC Fast-Block Correlation Report

## Goal

Check whether `fast_blocks=1` stays extremely close to dense sampling for:

- `BayesC`
- annotated `BayesC`

under same-seed paired production runs.

## Protocol

Script:

- `benchmarks/bayesc_fast_block_correlation.jl`

Settings:

- data seed: `20260401`
- MCMC seeds: `2026,2027`
- `n_obs = 50`
- `n_markers = 100`
- `chain_length = 800`
- `burnin = 200`
- `output_samples_frequency = 100`

Output directory:

- `/tmp/bayesc_fast_block_correlation_run`

Same-seed pairs:

- `BayesC_dense` vs `BayesC_fast_blocks_1`
- `Annotated_BayesC_dense` vs `Annotated_BayesC_fast_blocks_1`

## Pairwise Results

| Method | Seed | Marker Estimate Corr | Model Frequency Corr | EBV Corr | Residual Var Abs Diff | `pi_j` Corr | Annotation Coeff Corr |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `BayesC` | `2026` | `1.0` | `1.0` | `1.0` | `6.36e-7` | `1.0` | `NA` |
| `BayesC` | `2027` | `1.0` | `1.0` | `1.0` | `2.35e-6` | `1.0` | `NA` |
| `Annotated_BayesC` | `2026` | `1.0` | `1.0` | `1.0` | `1.29e-7` | `1.0` | `1.0` |
| `Annotated_BayesC` | `2027` | `1.0` | `1.0` | `1.0` | `4.57e-7` | `1.0` | `1.0` |

## Exact-Difference Check for Seed `2026`

Additional direct file-level checks on the saved outputs gave:

### BayesC

- max marker-effect absolute difference: `4.5e-7`
- max model-frequency absolute difference: `0.0`
- max EBV absolute difference: `4.49e-6`
- max `pi` output absolute difference: `0.0`

### Annotated BayesC

- max marker-effect absolute difference: `9.91e-7`
- max model-frequency absolute difference: `0.0`
- max EBV absolute difference: `5.14e-6`
- max `pi` output absolute difference: `0.0`
- max annotation-coefficient absolute difference: `0.0`

## Interpretation

For this current-master benchmark, both comparisons are effectively identical.

That is:

- standard `BayesC` with `fast_blocks=1` is not just highly correlated with
  dense `BayesC`; it is numerically almost identical here
- annotated `BayesC` with `fast_blocks=1` is also numerically almost identical
  to dense annotated `BayesC` here

So the current evidence from this same-seed production benchmark is:

- the user expectation is correct for these settings
- the previously suspected annotated-BayesC-specific instability is **not**
  reproduced on the current merged `master`

## Conclusion

On the current production code:

1. `BayesC` dense vs `fast_blocks=1` shows effectively exact same-seed agreement
   on this benchmark.
2. annotated `BayesC` dense vs `fast_blocks=1` also shows effectively exact
   same-seed agreement on this benchmark.
3. If a discrepancy still appears in a larger real analysis, it likely depends
   on a more extreme regime than the one reproduced here and should be debugged
   as a separate case rather than treated as a general property of annotated
   BayesC.
