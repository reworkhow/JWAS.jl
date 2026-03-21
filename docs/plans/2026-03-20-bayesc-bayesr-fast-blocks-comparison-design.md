# BayesC and BayesR Fast-Blocks Comparison Design

## Goal

Compare the production dense and `fast_blocks` paths for single-trait BayesC and
BayesR under one benchmark protocol.

## Scope

Benchmark matrix:

- BayesC dense
- BayesC `fast_blocks=true`
- BayesC `fast_blocks=1`
- BayesR dense
- BayesR `fast_blocks=true`
- BayesR `fast_blocks=1`

Protocol:

- single trait only
- production `runMCMC` path only
- long-chain, multiseed posterior-summary comparison
- no benchmark-local sampler loop as the main oracle

## Comparison Rule

Each method is compared against its own dense oracle:

- BayesC `fast_blocks=true` vs BayesC dense
- BayesC `fast_blocks=1` vs BayesC dense
- BayesR `fast_blocks=true` vs BayesR dense
- BayesR `fast_blocks=1` vs BayesR dense

The report should answer:

- how invariant BayesC is to block mode
- how invariant BayesR is to block mode
- whether `fast_blocks=1` is effectively identical to dense for each method
- how much drift default `fast_blocks=true` introduces for each method

## Implementation Shape

Keep this benchmark-only.

Add one benchmark script that:

- builds a shared synthetic single-trait dataset
- runs the six production cases above through `get_genotypes`, `build_model`,
  and `runMCMC`
- extracts method-specific posterior summaries
- writes:
  - per-run CSV
  - pairwise summary CSV
  - markdown report

Use a shared subset of metrics across methods:

- residual variance
- model or nonzero frequency
- runtime

Then keep method-specific metrics separate:

- BayesC: scalar `pi`, marker variance
- BayesR: 4-class `pi`, shared `sigmaSq`
