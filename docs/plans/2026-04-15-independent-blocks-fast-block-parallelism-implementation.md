# Independent-Blocks Fast-Block Parallelism Implementation

**Date:** 2026-04-15

**Branch:** `feature/independent-blocks-fast-blocks`

## Summary

This implementation adds an explicit `independent_blocks=true` mode for fast-block marker samplers.
The default remains `independent_blocks=false`, which keeps the existing exact sequential block sweep.

The new mode is intended for users who provide marker blocks that they are willing to treat as independent enough for parallel block updates.
Within each MCMC marker sweep, JWAS freezes the corrected phenotype / residual state, updates blocks independently, then applies all block effect deltas after the block barrier.

## Public API

`runMCMC` now accepts:

```julia
independent_blocks::Bool = false
```

Rules:

- `independent_blocks=false` keeps the exact sequential fast-block behavior.
- `independent_blocks=true` requires `fast_blocks != false`.
- `fast_blocks` now accepts explicit marker block starts, for example `fast_blocks=[1, 501, 975]`.
- Explicit block starts use full-sweep chain semantics: one MCMC iteration means one full sweep over the supplied blocks.
- `fast_blocks=true` and numeric `fast_blocks` keep the legacy fixed-size block behavior and chain-length scaling.

## Implemented Sampler Coverage

The independent-block path is implemented for:

- single-trait BayesA/B/C through `BayesABC_block!`
- single-trait BayesR through `BayesR_block!`
- single-trait annotated BayesC
- single-trait annotated BayesR
- plain multi-trait BayesA/B/C through `MTBayesABC_block!`
- annotated 2-trait BayesC

For multi-trait BayesC, the block sampler supports:

- `multi_trait_sampler=:I`
- `multi_trait_sampler=:II`
- `multi_trait_sampler=:auto`

## Algorithm

The exact block path is unchanged:

1. compute the block RHS from the current corrected phenotype / residual;
2. update markers inside the block;
3. update the global corrected phenotype / residual;
4. move to the next block.

The independent-block path changes only the inter-block schedule:

1. copy the sweep-level corrected phenotype / residual state;
2. generate block-local RNG seeds from the chain RNG;
3. update each block from the frozen state using block-local buffers and RNGs;
4. collect each block's old-minus-new marker-effect deltas;
5. after the threaded loop finishes, apply all block deltas to the global corrected phenotype / residual.

This avoids shared writes during block updates.
Julia `Threads.@threads` provides the block-level parallel execution.

## Statistical Assumption

For single-trait models, independent blocks are exact only when:

```text
X_b' W X_c = 0  for all b != c
```

For multi-trait models, trait coupling through `R^{-1}` remains inside the within-block sampler.
The independent-block approximation is still governed by off-block genotype leakage.

If off-block weighted crossproducts are not zero, `independent_blocks=true` is an explicit approximation for speed and parallelism.

## Tests Added

The test suite now covers:

- `independent_blocks` API acceptance.
- Clear error when `independent_blocks=true` is used with `fast_blocks=false`.
- Explicit block-start validation: empty, missing start at `1`, unsorted, duplicated, and out-of-range starts.
- Chain-length behavior for explicit block starts.
- Production-path single-trait BayesC independent blocks.
- Production-path annotated BayesC independent blocks.
- Production-path BayesR independent blocks.
- Production-path annotated BayesR independent blocks.
- Production-path plain multi-trait BayesC sampler I independent blocks.
- Production-path plain multi-trait BayesC sampler II independent blocks.
- Production-path annotated multi-trait BayesC sampler I independent blocks.
- Production-path annotated multi-trait BayesC sampler II independent blocks.
- A deterministic weighted crossproduct check documenting the condition under which independent blocks are exact.
- Two-thread execution smoke coverage for annotated and multi-trait independent-block paths.

Fresh verification passed:

- focused independent-block and annotated/multi-trait tests, exit code 0;
- two-thread focused independent-block tests, exit code 0;
- BayesR wrapper/probe regression test, exit code 0;
- docs build, exit code 0;
- full `test/runtests.jl`, `836/836` tests passing in `3m11.2s`.

## Benchmarks

The production smoke benchmark driver is:

```text
benchmarks/independent_blocks_fast_blocks_smoke.jl
```

It runs production `runMCMC` for:

- BayesC
- BayesR
- annotated BayesC
- annotated BayesR
- plain multi-trait BayesC sampler I
- plain multi-trait BayesC sampler II
- annotated multi-trait BayesC sampler I
- annotated multi-trait BayesC sampler II

It compares:

- exact sequential fast blocks: `independent_blocks=false`
- approximate independent blocks: `independent_blocks=true`

The current smoke outputs are saved under:

```text
benchmarks/reports/2026-04-15-independent-blocks-fast-block-parallelism-smoke-results-threads1.csv
benchmarks/reports/2026-04-15-independent-blocks-fast-block-parallelism-smoke-results-threads2.csv
```

These smoke runs are designed to verify production-path execution and record rough timing.
They are not a final performance benchmark because the demo dataset has only five markers and first-call compilation dominates some timings.

## Operational Guidance

Recommended server setup:

```bash
export JULIA_NUM_THREADS=<num_cores>
export OPENBLAS_NUM_THREADS=1
julia --project=.
```

Use `OPENBLAS_NUM_THREADS=1` to avoid nested thread oversubscription while JWAS parallelizes across marker blocks.

## Remaining Follow-Up

Future benchmarking should use a larger livestock genotype fixture and compare prediction accuracy, EBV correlations, marker-effect summaries, and posterior inclusion summaries.
Any marker-level comparison must align markers by explicit marker ID or SNP index before computing correlations.
