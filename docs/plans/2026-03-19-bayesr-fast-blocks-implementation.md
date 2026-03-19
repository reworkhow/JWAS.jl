# BayesR Fast Blocks Implementation Record

## Feature

Single-trait dense `fast_blocks` support for BayesR.

## Status

Partially implemented, not acceptance-ready.

The code path exists and runs, but the dense-vs-block posterior-summary benchmark
does not yet meet the design target for parity against dense BayesR.

## What Was Implemented

### Validation

Updated `input_data_validation.jl` so BayesR is allowed to use `fast_blocks`
while still rejecting the current unsupported combinations:

- `storage=:stream`
- multi-trait
- annotations
- RRM

### Sampler Integration

Added BayesR `fast_blocks` support in:

- `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

Implementation shape:

- dense BayesR still uses `BayesR!`
- BayesR with `fast_blocks != false` dispatches to `BayesR_block!`
- the block kernel mirrors the BayesABC block structure:
  - block RHS cache
  - per-marker updates within block
  - block-level `yCorr` correction on exit

### Tests

Extended `test/unit/test_bayesr.jl` with:

- validation coverage for BayesR `fast_blocks`
- unsupported-combination coverage
- direct block-kernel coverage
- dispatch coverage confirming `runMCMC(...; fast_blocks=true)` enters the
  BayesR block path

## Verification Run During Implementation

Targeted BayesR test command:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Observed result at the time of implementation:

- BayesR validation: pass
- BayesR dense sampler: pass
- BayesR block sampler: pass
- BayesR fast_blocks dispatch: pass

## Benchmark Harness Added

Added:

- `benchmarks/bayesr_fast_blocks_parity.jl`

Purpose:

- compare dense BayesR vs block BayesR inside JWAS
- report posterior summaries and runtime
- support an explicit numeric block size for diagnostics

Important harness note:

- the benchmark uses a block-aware burnin for the block run because JWAS rescales
  `chain_length` internally under `fast_blocks`

## What Was Discovered

The first parity run initially showed all-zero block summaries.

Root cause:

- harness bug, not sampler bug
- `burnin` exceeded the internally rescaled block outer chain length, so no
  post-burnin samples were written

After fixing the harness, the real issue became clear:

- BayesR block runs are substantially faster
- but posterior summaries do not yet match dense BayesR closely enough

Key observations from the saved parity report:

- default block size from `fast_blocks=true` is too aggressive in the benchmarked
  setting
- smaller explicit block sizes improve parity
- even with smaller blocks, the current BayesR block path still shows material
  drift in `sigmaSq`, `pi`, and nonzero frequency

See:

- `benchmarks/reports/2026-03-19-bayesr-fast-blocks-parity.md`

## Local Hypothesis Checks Performed

One local schedule experiment was tried and discarded:

- changing BayesR block `nreps` from `block_size` to `1`

Result:

- did not improve parity enough
- was reverted

So the current mismatch does not look like a one-line schedule-constant issue.

## Current Recommendation

Do not treat BayesR `fast_blocks` as merge-ready yet.

Recommended next step:

1. benchmark dense vs block with `estimatePi=false`
2. isolate whether the main drift enters through:
   - block marker updates
   - class occupancy
   - `sigmaSq` update
3. reconsider the BayesR block schedule rather than assuming BayesABC block
   scheduling carries over unchanged

## Files Touched In This Iteration

- `src/1.JWAS/src/input_data_validation.jl`
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- `test/unit/test_bayesr.jl`
- `benchmarks/bayesr_fast_blocks_parity.jl`
- `benchmarks/reports/2026-03-19-bayesr-fast-blocks-parity.md`
