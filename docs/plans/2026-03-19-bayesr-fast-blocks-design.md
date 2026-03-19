# BayesR Fast Blocks Design

**Date:** 2026-03-19

## Goal

Add a first `fast_blocks` implementation for BayesR in JWAS, scoped to single-trait dense genotype analysis and validated against the existing dense BayesR path by posterior-summary parity.

## Scope

### Implemented in this design
- Single-trait `method="BayesR"` support for `fast_blocks=true` or numeric block size.
- Dense genotype storage only.
- Reuse of the existing BayesR model specification:
  - fixed `BAYESR_GAMMA = [0.0, 0.01, 0.1, 1.0]`
  - shared marker variance `sigmaSq`
  - 4-class mixture probabilities `pi`
- Validation and tests against the current dense BayesR implementation.

### Explicitly out of scope
- `storage=:stream`
- Multi-trait BayesR
- BayesR annotations
- RRM
- New public API specific to BayesR block updates

## Architecture

BayesR `fast_blocks` will follow the same structural split already used by BayesC:
- keep the current dense BayesR kernel in `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- add a new `BayesR_block!` wrapper and a lower-level block kernel in the same file
- dispatch to the block kernel from `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl` when:
  - `Mi.method == "BayesR"`
  - single-trait
  - dense storage
  - `fast_blocks != false`

The outer MCMC flow should remain unchanged:
- same setup of block matrices through `GibbsMats(..., fast_blocks=...)`
- same `pi` update
- same marker-variance update
- same output semantics

Validation changes should be minimal:
- remove the current BayesR rejection of `fast_blocks`
- continue rejecting `storage=:stream`, multi-trait, annotations, and RRM for BayesR

## Block Sampler Logic

The BayesR block kernel should mirror `BayesABC_block!` structurally while keeping BayesR-specific class sampling.

Per block:

1. Build the cached block right-hand side once using `block_rhs!`.
2. Save the current marker effects for the block in `αold_block`.
3. Within the block inner loop, update one marker at a time:
   - compute the BayesR marker `rhs` from the cached block state
   - evaluate the 4 class log probabilities
   - sample the class label
   - sample a new marker effect if the class is nonzero
   - update the cached block RHS with `oldAlpha - newAlpha`
4. After the block finishes, apply one block-level correction back to `yCorr` using `mul!`, as in the BayesABC block path.

The first implementation should keep the existing BayesABC block schedule:
- `nreps = block_size`

The design target is code-structure parity with BayesABC, not a new BayesR-specific block schedule.

## Validation Strategy

Validation should treat dense BayesR as the oracle.

### Required checks
- class labels stay in `1:4`
- `pi` remains a valid probability vector
- unsupported combinations still error:
  - `storage=:stream`
  - multi-trait
  - annotations
  - RRM

### Parity target

Compare dense BayesR vs `fast_blocks` BayesR by posterior summaries, not iteration-level identity:
- posterior mean `sigmaSq`
- posterior mean residual variance
- posterior mean `pi`
- posterior mean nonzero frequency
- posterior mean marker-effect correlation

The acceptance target is posterior-summary parity, not exact trajectory matching.

### Performance target

Record at least one runtime comparison on a moderate dense dataset. Runtime improvement is useful, but parity against dense BayesR is the primary acceptance criterion.

## Files To Modify

- `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
  - add `BayesR_block!`
  - add any small BayesR block helper(s)
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
  - dispatch BayesR to the block path when `fast_blocks != false`
- `src/1.JWAS/src/input_data_validation.jl`
  - allow BayesR `fast_blocks`
  - keep other BayesR restrictions
- `test/unit/test_bayesr.jl`
  - add validation and output tests for `fast_blocks`
- benchmark and report files as needed for dense-vs-block parity

## Deliverables

- Code support for single-trait dense BayesR `fast_blocks`
- Targeted tests
- Full JWAS test-suite pass
- A benchmark/report artifact showing posterior-summary parity between dense and block BayesR
