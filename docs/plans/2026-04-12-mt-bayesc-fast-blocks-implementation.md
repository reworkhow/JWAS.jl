# Fast-Block Multi-Trait BayesC Implementation

Date: 2026-04-12

## Scope Delivered

- Changed the default `multi_trait_sampler` for multi-trait BayesC from `:auto` to `:I`, while keeping explicit `:auto` and `:II`.
- Refactored the multi-trait BayesC block wrapper so dense and fast-block paths both choose:
  - prior source via `GlobalPiPrior` or `MarkerSpecificPiPrior`
  - sampler mode via `mt_bayesc_sampler_mode(...)`
- Implemented fast-block sampler I and sampler II for multi-trait BayesC.
- Enabled annotated 2-trait BayesC with `fast_blocks=true`.
- Removed the old `fast_blocks=false` restriction for annotated multi-trait BayesC, while keeping the existing restrictions on:
  - exactly 2 traits
  - `storage=:dense`
  - `constraint=false`
- Hardened sampler-II probability normalization so explicit `:II` runs remain valid even when some states have zero prior mass.

## Key Code Changes

- `src/1.JWAS/src/markers/readgenotypes.jl`
  - default `multi_trait_sampler` is now `:I`
- `src/1.JWAS/src/types.jl`
  - genotype default updated accordingly
- `src/1.JWAS/src/build_MME.jl`
  - explicit unsupported override checks remain targeted to `:II`
- `src/1.JWAS/src/input_data_validation.jl`
  - removed the old `fast_blocks=false` block for annotated multi-trait BayesC
  - removed the old `fast_blocks=false` block for explicit sampler overrides
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
  - runtime guard no longer rejects annotated multi-trait BayesC in block mode
- `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
  - added `mt_bayesc_block_context(...)`
  - block wrapper now dispatches by prior source and sampler mode
  - added `_MTBayesABC_block_samplerI!`
  - added `_MTBayesABC_block_samplerII!`
  - sampler-II probability normalization now uses a stable softmax-style normalization in both dense and block paths

## Test Coverage Added or Updated

- `test/unit/test_multitrait_mcmc.jl`
  - default multi-trait BayesC sampler is `:I`
  - explicit `multi_trait_sampler=:auto` on the plain multi-trait path still follows the legacy support-based rule
  - explicit dense sampler-II run succeeds for plain multi-trait BayesC
  - fast-block sampler-I run succeeds for plain multi-trait BayesC
  - fast-block sampler-II run succeeds for plain multi-trait BayesC
  - fast-block sampler-II run succeeds for annotated multi-trait BayesC
  - block wrapper chooses the expected prior source and sampler mode
  - one-marker block sampler II matches dense sampler II on the same target posterior
- `test/unit/test_annotated_bayesc.jl`
  - default multi-trait annotated BayesC sampler is `:I`
  - explicit `:auto` is preserved
  - annotated multi-trait BayesC now runs with `fast_blocks=true`

## Verification

Focused:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Full suite:

```bash
julia --project=. --startup-file=no test/runtests.jl
```
