# BayesC And BayesR Fast-Blocks EBV Correlation Design

## Goal

Extend the production fast-block comparison benchmark so it also reports in-sample
prediction accuracy as the correlation between observed phenotype `y` and posterior
mean EBV.

## Scope

- single-trait only
- production `runMCMC` path only
- benchmark script:
  - `benchmarks/bayesc_bayesr_fast_blocks_comparison.jl`
- methods already in the benchmark:
  - `BayesC_dense`
  - `BayesC_fast_blocks_default`
  - `BayesC_fast_blocks_1`
  - `BayesR_dense`
  - `BayesR_fast_blocks_default`
  - `BayesR_fast_blocks_1`

## Metric

For each method variant and seed:

1. request EBV output from production `runMCMC`
2. read `output["EBV_y1"]`
3. join EBV with the synthetic phenotype table on `ID`
4. compute:
   - `cor(y1, EBV)`

This metric will be stored as `phenotype_ebv_correlation`.

## Output Changes

- add `phenotype_ebv_correlation` to `comparison_runs.csv`
- add pairwise dense-vs-block absolute differences for this metric to
  `comparison_pairwise_summary.csv`

## Testing

Update the benchmark subprocess test so it fails until:

- `comparison_runs.csv` contains `phenotype_ebv_correlation`
- `comparison_pairwise_summary.csv` contains metric
  `phenotype_ebv_correlation`
