# BayesR Fast Blocks Burnin Repetition Schedule Implementation

**Date:** 2026-03-19

## Summary

Implemented a BayesR-specific fast-block scheduling rule:

- during effective burnin, BayesR block updates use `nreps = 1`
- after effective burnin, BayesR block updates use `nreps = block_size`

This changes only BayesR fast blocks. Dense BayesR and BayesABC block behavior are unchanged.

## Code Changes

### Production

- [`src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`](/Users/haocheng/Github/JWAS.jl/src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl)
  - added `bayesr_block_nreps(iter, burnin, block_size)`
  - added `BayesR_block!` methods that accept `iter` and `burnin`
  - replaced the hardcoded block repetition rule with the helper

- [`src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`](/Users/haocheng/Github/JWAS.jl/src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl)
  - BayesR fast-block calls now pass the current `iter` and effective `burnin` into `BayesR_block!`

### Tests

- [`test/unit/test_bayesr.jl`](/Users/haocheng/Github/JWAS.jl/test/unit/test_bayesr.jl)
  - added a regression test for `bayesr_block_nreps`
  - updated the BayesR block-dispatch probe to match the new call signature

## TDD Record

Red:

- added the `BayesR fast-block repetition schedule` testset
- ran:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

- failure was the expected one:
  - `UndefVarError: bayesr_block_nreps not defined in JWAS`

Green:

- implemented the helper and threaded `iter` / `burnin` into the BayesR block path
- reran the same targeted test file and it passed

## Verification

Targeted:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Result:

- all BayesR testsets passed, including the new scheduling regression

Full suite:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Result:

- `543/543` passing

## Focused Behavior Check

I also ran one focused production comparison against dense BayesR using:

- same dataset
- same seed
- free `pi`
- free `sigmaSq`
- block-aware burnin

Results after the burnin-gated change:

- default `fast_blocks=true`
  - `sigmaSq` abs diff: `0.0580`
  - residual variance abs diff: `0.0329`
  - mean nonzero frequency abs diff: `0.0538`

- `fast_blocks=2`
  - `sigmaSq` abs diff: `0.0935`
  - residual variance abs diff: `0.0103`
  - mean nonzero frequency abs diff: `0.0176`

The main occupancy / `pi` mismatch is much smaller than before for both default blocks and block size 2. The tradeoff is that `sigmaSq` did not uniformly improve in the same direction, especially for `block_size=2`.

So this implementation clearly changes the BayesR block dynamics in the intended direction for early occupancy stabilization, but it does not close every summary-level gap.

## Local-Only Files Left Out Of This Change

These earlier exploratory benchmark files remain local-only and are intentionally not part of the narrow production change:

- [`benchmarks/bayesr_fast_blocks_parity.jl`](/Users/haocheng/Github/JWAS.jl/benchmarks/bayesr_fast_blocks_parity.jl)
- [`benchmarks/bayesr_parity_common.jl`](/Users/haocheng/Github/JWAS.jl/benchmarks/bayesr_parity_common.jl)
- [`test/unit/test_bayesr_parity.jl`](/Users/haocheng/Github/JWAS.jl/test/unit/test_bayesr_parity.jl)
