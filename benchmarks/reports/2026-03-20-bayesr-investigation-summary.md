# BayesR Investigation Summary

## Executive Summary

The BayesR investigation has now answered the main questions clearly.

### Dense BayesR vs the R reference

- There is no evidence of a substantive BayesR logic bug in dense JWAS.
- One-step controlled replay matched the R reference under shared draws.
- `Float32` precision was far too small an effect to explain the original gap.
- Production JWAS and benchmark-local Julia BayesR agreed closely.
- The earlier dense JWAS vs R mismatch was mostly short-chain Monte Carlo noise.
- Under the final long-chain multiseed protocol, JWAS and the R reference agree
  well at the posterior-summary level.

### BayesR `fast_blocks`

- The block machinery itself is sound.
- `fast_blocks=1` matches dense BayesR essentially exactly with the same seed.
- Default blocks with `nreps=1` also match dense BayesR essentially exactly in
  the benchmark-local comparison.
- The remaining drift enters when repeated within-block sweeps are used.
- A burnin-gated schedule improves the behavior materially.
- On long chains, production `fast_blocks` is much closer to dense BayesR than
  the early short-chain tests suggested.
- A smaller explicit block size (`2`) improves parity further, at a speed cost.

## Part 1: Dense BayesR Parity Against R

### 1. Controlled replay ruled out a one-step formula bug

From [2026-03-18-bayesr-controlled-replay-report.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-18-bayesr-controlled-replay-report.md):

- same exported draws were fed to JWAS-local replay and the R reference
- `chosen_class` matched for every marker
- `new_alpha` matched for every marker
- `ssq` matched exactly
- `nnz` matched exactly
- `sigmaSq_new` matched to numerical precision (`2.27e-12`)

Meaning:

- the BayesR update formulas align with the R reference
- the original iteration-1 divergence was not evidence of an immediate logic bug

### 2. Float precision was not the explanation

From [2026-03-18-bayesr-parity-monte-carlo-report.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-18-bayesr-parity-monte-carlo-report.md):

- `Float32` vs `Float64` inside Julia, with shared draws:
  - `sigmaSq` relative diff: `3.47e-8`
  - residual-variance relative diff: `3.43e-7`

Meaning:

- the short-chain parity gap was not caused by `Float32`

### 3. Production JWAS and local Julia BayesR agreed

Also from [2026-03-18-bayesr-parity-monte-carlo-report.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-18-bayesr-parity-monte-carlo-report.md):

- same dataset
- same seed
- fixed `pi`
- `chain_length=1000`
- `burnin=200`

Result:

- production JWAS `sigmaSq`: `6.329144`
- benchmark-local Julia `sigmaSq`: `6.327308`
- relative difference: about `0.03%`

Meaning:

- the remaining parity question was not caused by a production-vs-local Julia gap

### 4. Short-chain dense parity was misleading

From [2026-03-18-bayesr-parity-monte-carlo-report.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-18-bayesr-parity-monte-carlo-report.md):

For fixed `pi`, `1000/200`, across seeds `2021:2030`:

- mean JWAS-vs-R `sigmaSq` relative difference: `9.93%`
- max `sigmaSq` relative difference: `27.40%`

But within-method spread was already large:

- JWAS `sigmaSq` mean `7.8340`, SD `2.0567`
- R `sigmaSq` mean `8.7269`, SD `2.8776`

Meaning:

- short single-chain parity was too noisy to use as a pass/fail rule

### 5. Long-chain multiseed dense parity was good

From [2026-03-18-bayesr-parity-final-note.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-18-bayesr-parity-final-note.md):

Protocol:

- `chain_length=10000`
- `burnin=2000`
- seeds `2026:2030`

Fixed `pi`:

- mean `sigmaSq` relative difference: `0.032145`
- max `sigmaSq` relative difference: `0.055008`
- mean residual-variance relative difference: `0.005274`
- mean nonzero-frequency absolute difference: `0.001446`

`estimate_pi`:

- mean `sigmaSq` relative difference: `0.015964`
- max `sigmaSq` relative difference: `0.043060`
- mean residual-variance relative difference: `0.006239`
- mean nonzero-frequency absolute difference: `0.007581`
- mean max-`pi` absolute difference: `0.011302`

Meaning:

- dense JWAS BayesR is in good shape relative to the R reference

## Part 2: BayesR `fast_blocks` Investigation

### 1. Initial short-chain block results looked bad

From [2026-03-19-bayesr-fast-blocks-parity.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-19-bayesr-fast-blocks-parity.md):

With early short-chain runs:

- default block size `7` showed large drift in `sigmaSq`, occupancy, and `pi`
- block size `2` improved the result, but still showed material drift

That result was real, but it was not the full story.

### 2. Two real benchmark bugs were found and fixed

#### `fast_blocks=1` was being mistaken for `true`

Cause:

- code used `fast_blocks == true`
- in Julia, `1 == true` is `true`

Effect:

- `fast_blocks=1` silently resolved to the default block size instead of block size `1`

Fix:

- committed in `9d4aa267`
- numeric `fast_blocks` values are now distinguished from booleans correctly

#### `Pi` input aliasing caused same-session benchmark contamination

Cause:

- `get_genotypes` stored the caller’s `Pi` vector by reference
- BayesR updates `Mi.π` in place
- a first run mutated the starting `Pi` for the second run

Effect:

- a supposed same-seed dense-vs-block comparison was not actually starting from
  the same state

Fix:

- committed in `72ceb48f`
- `get_genotypes` now copies mutable `Pi` input vectors

### 3. `fast_blocks=1` matches dense BayesR

After fixing the two benchmark bugs:

- same seed
- free `pi`
- free `sigmaSq`
- `fast_blocks=1`

Result:

- dense and block BayesR matched to numerical precision

Meaning:

- the basic block kernel is not wrong

### 4. Default blocks with `nreps=1` also match dense in the benchmark-local comparison

This was tested with:

- default block partition
- one sweep per SNP inside each block
- same outer chain length as dense

Result:

- dense vs `default_blocks_single_rep` matched essentially exactly on:
  - `sigmaSq`
  - residual variance
  - nonzero frequency
  - `pi`
  - marker effects

Meaning:

- the block partition itself is not the main problem
- the repeated within-block sweeps are the main source of drift

### 5. Fixed-hyperparameter diagnostics narrowed the problem

From [2026-03-19-bayesr-fast-blocks-fixed-hyperparameters.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-19-bayesr-fast-blocks-fixed-hyperparameters.md):

With:

- fixed `pi`
- fixed `sigmaSq`
- residual variance still sampled

Result:

- block-vs-dense nonzero-frequency abs diff: `0.000339`
- dense-only seed range for nonzero frequency: `0.006500`
- residual-variance abs diff: `0.038778`
- dense-only residual-variance seed range: `0.013856`

Meaning:

- most of the earlier occupancy drift came from the `pi` / `sigmaSq` feedback loop
- the remaining mismatch was much narrower than the original short-chain result

### 6. Production burnin-gated long-chain BayesR fast blocks look much better

Fresh long-chain production results from [2026-03-20-bayesr-fast-blocks-long-chain-schedule-report.md](/Users/haocheng/Github/JWAS.jl/benchmarks/reports/2026-03-20-bayesr-fast-blocks-long-chain-schedule-report.md):

Protocol:

- `chain_length=10000`
- `burnin=2000`
- seeds `2026:2030`
- free `pi`
- free `sigmaSq`

#### Default block size (`7`)

Burnin-gated vs dense:

- mean `sigmaSq` abs diff: `0.043457`
- mean residual-variance abs diff: `0.004309`
- mean nonzero-frequency abs diff: `0.036223`
- mean `pi_class1` abs diff: `0.033362`
- speedup vs dense: `15.29x`

#### Explicit block size (`2`)

Burnin-gated vs dense:

- mean `sigmaSq` abs diff: `0.028267`
- mean residual-variance abs diff: `0.005784`
- mean nonzero-frequency abs diff: `0.022867`
- mean `pi_class1` abs diff: `0.021246`
- speedup vs dense: `9.16x`

Meaning:

- the long-chain production behavior is much better than the early short-chain
  diagnostics suggested
- a smaller block size improves most BayesR-specific summaries
- there is a clear speed-vs-parity tradeoff

## What Is A Bug And What Is Not

### Confirmed bugs

- `fast_blocks=1` treated as boolean `true`
- `Pi` aliasing in `get_genotypes`

### Not supported as bugs by the evidence

- dense BayesR logic bug
- BayesR one-step formula mismatch vs R
- `Float32` precision as the explanation for the original dense parity gap
- basic BayesR block-kernel algebra bug

### What remains

The remaining fast-block question is not basic correctness. It is algorithmic and
practical:

- how much repeated within-block sweeping should BayesR use
- and what block size is acceptable for the desired speed/parity tradeoff

## Current Recommendations

### Dense BayesR

- treat dense BayesR as validated
- use long-chain or multiseed protocols for future parity work
- do not use short single chains as decisive evidence

### BayesR `fast_blocks`

- do not judge BayesR fast blocks from short-chain same-seed traces
- if parity is the main goal:
  - prefer smaller block sizes
  - or keep BayesR at one sweep per SNP
- if speed is the main goal:
  - the burnin-gated production schedule is defensible
  - but it is an approximation tradeoff, not the same kernel as dense BayesR

## Bottom Line

The investigation changed the story substantially.

The initial suspicion was:

- dense BayesR might be wrong
- BayesR `fast_blocks` might be fundamentally broken

The evidence now supports a different conclusion:

- dense BayesR in JWAS is in good shape
- the original dense parity gap was mostly a benchmarking/protocol issue
- BayesR `fast_blocks` is not fundamentally broken
- the real remaining choice is how aggressively BayesR should exploit repeated
  within-block sweeps, and what speed/parity tradeoff is acceptable
