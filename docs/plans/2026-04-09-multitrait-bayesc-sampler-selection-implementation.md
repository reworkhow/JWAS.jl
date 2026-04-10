# Multi-Trait BayesC Sampler Selection Implementation

**Date:** 2026-04-09

## Summary

Added an explicit `multi_trait_sampler` option to `get_genotypes` and `Genotypes` so multi-trait BayesC can choose:

- `:auto` — preserve the existing support-based dispatch
- `:I` — force Gibbs sampler I
- `:II` — force Gibbs sampler II

The default remains `:auto`, so existing workflows keep the same behavior unless the user opts in.

## Public API

`get_genotypes(...; multi_trait_sampler=:auto)` now accepts:

- `:auto`
- `:I`
- `:II`

Invalid symbols fail immediately in `get_genotypes`.

## Dispatch Behavior

`MTBayesABC!` now uses `mt_bayesc_sampler_mode(genotypes, nModels)` to decide between sampler I and II:

- `:auto` keeps the legacy rule `length(genotypes.π) == 2^nModels ? :I : :II`
- `:I` forces `_MTBayesABC_samplerI!`
- `:II` forces `_MTBayesABC_samplerII!`

This applies to both:

- ordinary multi-trait BayesC with a global `Pi` dict
- annotated 2-trait BayesC with marker-specific `snp_pi`

The prior-source logic (`GlobalPiPrior` vs `MarkerSpecificPiPrior`) is unchanged.

## Validation

Explicit sampler overrides are now validated in two stages:

1. `get_genotypes`
   - rejects invalid symbols

2. `build_model` / runtime validation
   - rejects explicit overrides outside multi-trait BayesC
   - keeps existing annotated multi-trait BayesC restrictions unchanged
   - still requires `storage=:dense` and `fast_blocks=false` for annotated multi-trait BayesC

This means annotated 2-trait BayesC can now intentionally use sampler II, but only through the already supported dense non-block path.

## Tests

Added/updated coverage for:

- invalid `multi_trait_sampler` symbols
- storing `multi_trait_sampler` on annotated BayesC genotype objects
- rejecting explicit overrides in single-trait BayesC builds
- forcing sampler II in annotated 2-trait BayesC and running the production MCMC path end to end
