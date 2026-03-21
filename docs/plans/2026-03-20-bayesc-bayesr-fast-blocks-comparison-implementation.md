# BayesC and BayesR Fast-Blocks Comparison Implementation

## Summary

Implemented a benchmark-only production comparison for single-trait BayesC and
BayesR dense vs `fast_blocks`.

## Files

- added:
  [bayesc_bayesr_fast_blocks_comparison.jl](/Users/haocheng/Github/JWAS.jl/benchmarks/bayesc_bayesr_fast_blocks_comparison.jl)
- modified:
  [test_bayesr_parity.jl](/Users/haocheng/Github/JWAS.jl/test/unit/test_bayesr_parity.jl)
- added report:
  [2026-03-20-bayesc-bayesr-fast-blocks-comparison-report.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-20-bayesc-bayesr-fast-blocks-comparison-report.md)

## Benchmark Protocol

The benchmark runs six production cases:

- BayesC dense
- BayesC `fast_blocks=true`
- BayesC `fast_blocks=1`
- BayesR dense
- BayesR `fast_blocks=true`
- BayesR `fast_blocks=1`

The script uses:

- a shared synthetic dataset
- multiseed runs
- block-aware burnin for the block paths
- production `get_genotypes`, `build_model`, and `runMCMC`

## Verification

Ran:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'`
- `julia --project=. --startup-file=no benchmarks/bayesc_bayesr_fast_blocks_comparison.jl /tmp/bayesc_bayesr_fast_blocks_comparison_20260320`
- `julia --project=. --startup-file=no test/runtests.jl`
- `julia --project=docs --startup-file=no docs/make.jl`

## Notes

- The raw runtime values in the benchmark output are affected by JIT warmup and
  should be interpreted cautiously.
- The benchmark was intentionally kept out of the production sampler path.
