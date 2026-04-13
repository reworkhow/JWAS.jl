# Fast-Block Multi-Trait BayesC Design

## Goal

Extend the multi-trait BayesC `fast_blocks` path so it supports the same sampler split and prior abstractions as the current dense path:

- plain multi-trait BayesC
- annotated 2-trait multi-trait BayesC
- sampler I
- sampler II
- explicit `multi_trait_sampler = :I | :II | :auto`

The implementation should preserve the current dense semantics and make the block path a block-form translation of the same samplers, rather than a separate algorithm.

## Current State

The dense multi-trait BayesC path in [MTBayesABC.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/mt-bayesc-fast-blocks/src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl) is already structured around:

- `GlobalPiPrior`
- `MarkerSpecificPiPrior`
- `log_prior(prior, marker, state)`
- `mt_bayesc_sampler_mode(...)`
- `_MTBayesABC_samplerI!`
- `_MTBayesABC_samplerII!`

The block path is still the older implementation:

- one sampler-I-style kernel only
- direct use of `genotypes.π`
- no marker-specific prior access
- no respect for `multi_trait_sampler`

That is why:

- plain multi-trait BayesC block mode does not yet share the dense sampler split
- annotated multi-trait BayesC block mode is explicitly rejected

## Design Decision

The block path should not invent new model logic. It should be the same sampler logic as the dense path, translated into block-local residual algebra.

The design is:

1. Keep the dense multi-trait BayesC path unchanged.
2. Replace the old `MTBayesABC_block!` body with a wrapper that:
   - chooses prior source using `GlobalPiPrior` or `MarkerSpecificPiPrior`
   - chooses sampler mode using the same `mt_bayesc_sampler_mode(...)` helper as dense
   - dispatches to block sampler I or block sampler II
3. Implement:
   - `_MTBayesABC_block_samplerI!`
   - `_MTBayesABC_block_samplerII!`
4. Reuse:
   - `GlobalPiPrior`
   - `MarkerSpecificPiPrior`
   - `log_prior(...)`

## Sampler Policy

`multi_trait_sampler` keeps three values:

- `:I`
- `:II`
- `:auto`

The default user-facing choice will become `:I` for multi-trait BayesC, dense and block.

`multi_trait_sampler = :auto` remains available explicitly and keeps the existing dense support-based rule:

- sampler I if `length(genotypes.π) == 2^nModels`
- sampler II otherwise

This same rule will be used in the new block wrapper.

## Fast-Block Scope

This single refactor should cover both:

- plain multi-trait BayesC fast blocks
- annotated 2-trait multi-trait BayesC fast blocks

Once the block path can consume `MarkerSpecificPiPrior`, the current annotated multi-trait `fast_blocks=false` restriction can be removed.

Other current restrictions stay as they are unless explicitly required by the block implementation:

- annotated multi-trait BayesC remains exactly 2 traits
- `storage=:dense`
- `constraint=false`

## Testing Strategy

Coverage should prove both dispatch and correctness boundaries:

1. Dense multi-trait BayesC still respects `:I`, `:II`, and explicit `:auto`.
2. Default multi-trait BayesC now resolves to sampler I.
3. Plain multi-trait BayesC fast blocks run under:
   - default `:I`
   - explicit `:II`
   - explicit `:auto`
4. Annotated multi-trait BayesC fast blocks run under:
   - default `:I`
   - explicit `:II`
   - explicit `:auto`
5. The old validation rejecting annotated multi-trait `fast_blocks=true` is removed and replaced with coverage for the supported path.

The exact tiny-case sampler-I vs sampler-II target-equivalence tests remain useful and should be kept unchanged unless the new block path requires additional block-specific exact checks.

## Non-Goals

This change does not attempt to:

- unify dense and block kernels into one implementation
- change single-trait BayesC or BayesR block logic
- relax multi-trait annotated BayesC beyond the current dense-only assumptions except for `fast_blocks`
- redesign the `:auto` rule itself

## Expected Outcome

After this feature:

- dense and block multi-trait BayesC share the same sampler vocabulary
- plain and annotated multi-trait BayesC share the same prior abstraction in both dense and block paths
- annotated multi-trait BayesC gains a real `fast_blocks` production path
- the block path becomes reviewable in the same conceptual frame as the dense path
