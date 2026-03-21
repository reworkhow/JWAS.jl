# BayesR Fast-Blocks Benchmark Cleanup Implementation

## Summary

Cleaned up the BayesR fast-blocks benchmark layer without changing the
production BayesR sampler.

## Changes

### Main benchmark script

Kept [benchmarks/bayesr_fast_blocks_parity.jl](/Users/haocheng/Github/JWAS.jl/benchmarks/bayesr_fast_blocks_parity.jl)
focused on production-facing modes only:

- `fixed_hyperparameters`
- `dense_multiseed`
- `long_chain_schedule_comparison`

Removed benchmark-local local-loop helpers and debug-only modes from this file.

Also wrapped the script entrypoint in a `main()` function guarded by
`abspath(PROGRAM_FILE) == @__FILE__` so the script can be safely included by the
debug script.

### Debug benchmark script

Created [benchmarks/debug/bayesr_fast_blocks_debug.jl](/Users/haocheng/Github/JWAS.jl/benchmarks/debug/bayesr_fast_blocks_debug.jl)
for exploratory local-loop diagnostics.

This script now owns:

- local BayesR state initialization
- fixed-hyperparameter trace mode
- default-blocks single-rep mode

### Tests

Updated [test_bayesr_parity.jl](/Users/haocheng/Github/JWAS.jl/test/unit/test_bayesr_parity.jl)
to reflect the split:

- the production-facing long-chain schedule benchmark test now expects only:
  - `dense`
  - `burnin_gated`
- the debug-only single-rep subprocess test now targets the new debug script
- the generic trace comparator test remains in place

## Result

The benchmark layer is cleaner:

- production-facing benchmark modes stay in the main script
- exploratory local-loop diagnostics live in a separate debug script
- production BayesR code is unchanged

## Verification

Ran:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'`
- `julia --project=. --startup-file=no test/runtests.jl`
- `julia --project=docs --startup-file=no docs/make.jl`
