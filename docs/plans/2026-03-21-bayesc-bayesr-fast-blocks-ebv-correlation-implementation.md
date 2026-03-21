# BayesC And BayesR Fast-Blocks EBV Correlation Implementation

## Summary

Extended the production fast-block comparison benchmark so it also reports
in-sample prediction accuracy as the correlation between phenotype `y1` and
posterior mean `EBV`.

## Files

- modified:
  [bayesc_bayesr_fast_blocks_comparison.jl](/Users/haocheng/Github/JWAS.jl/benchmarks/bayesc_bayesr_fast_blocks_comparison.jl)
- modified:
  [test_bayesr_parity.jl](/Users/haocheng/Github/JWAS.jl/test/unit/test_bayesr_parity.jl)
- updated report:
  [2026-03-20-bayesc-bayesr-fast-blocks-comparison-report.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-20-bayesc-bayesr-fast-blocks-comparison-report.md)

## Implementation

- requested production EBV output from `runMCMC`
- joined `output["EBV_y1"]` with the phenotype table on `ID`
- computed `cor(y1, EBV)` as `phenotype_ebv_correlation`
- added the metric to:
  - `comparison_runs.csv`
  - `comparison_pairwise_summary.csv`

## Verification

Ran:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'`
- `julia --project=. --startup-file=no benchmarks/bayesc_bayesr_fast_blocks_comparison.jl /tmp/bayesc_bayesr_fast_blocks_comparison_ebv_20260321`

## Notes

- this is in-sample correlation, not cross-validation accuracy
- the benchmark still uses the production `runMCMC` path only
