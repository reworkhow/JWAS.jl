# BayesR Fast Blocks Burnin Repetition Schedule Design

**Date:** 2026-03-19

## Goal

Change single-trait dense BayesR `fast_blocks` so early block updates behave like dense BayesR, then switch to the current repeated within-block sweep after burnin.

## Problem

Current BayesR `fast_blocks` uses:

- `nreps = block_size`

for every outer iteration in [`src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`](/Users/haocheng/Github/JWAS.jl/src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl).

Diagnostics established:

- `fast_blocks=1` matches dense BayesR essentially exactly, including free `pi` and free `sigmaSq`, when the runs start from the same state.
- divergence appears when `block_size > 1`
- the main shift is in class occupancy and `pi`, not in an obvious one-step formula bug

So the risky part is not the block cache machinery itself. The risky part is repeated within-block sweeps during early chain stabilization.

## Intended Behavior

For BayesR only:

- if `fast_blocks == false`, keep the dense BayesR path unchanged
- if `fast_blocks != false`, keep the current block partitioning
- during effective burnin:
  - use `nreps = 1`
- after effective burnin:
  - use `nreps = block_size`

This makes early fast-block BayesR equivalent in update count to dense BayesR while retaining the existing repeated block sweeps for retained samples.

## Transition Rule

Let:

- `outer_iter` be the current fast-block outer iteration
- `burnin_outer` be the effective burnin used by `runMCMC` after any fast-block rescaling

Then:

- if `outer_iter <= burnin_outer`, use `nreps = 1`
- otherwise use `nreps = block_size`

This keeps the rule simple and tied to the existing MCMC semantics. No new public user option is added in the first cut.

## Scope

In scope:

- single-trait BayesR
- dense backend
- `fast_blocks != false`
- production `runMCMC` path

Out of scope:

- BayesABC block behavior
- multi-trait BayesR
- `storage=:stream`
- annotations
- new public API knobs for block scheduling

## Implementation Shape

The change should stay small:

- [`src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`](/Users/haocheng/Github/JWAS.jl/src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl)
  - make `BayesR_block!` accept the current repetition count instead of always deriving `nreps = block_size`
- [`src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`](/Users/haocheng/Github/JWAS.jl/src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl)
  - compute the BayesR fast-block repetition count from the current outer iteration and effective burnin
  - pass that repetition count into `BayesR_block!`

No global flags and no extra scheduler object are needed.

## Validation

Required validation:

- targeted BayesR unit tests
- regression test for the new scheduling rule
- full `test/runtests.jl`

Important benchmark expectation:

- `fast_blocks=1` remains unchanged and stays dense-equivalent
- `block_size > 1` should move closer to dense BayesR in early-chain behavior because burnin no longer uses repeated within-block sweeps

