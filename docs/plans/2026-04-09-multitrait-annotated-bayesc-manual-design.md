# Multi-Trait Annotated BayesC Manual Design

## Goal

Document dense 2-trait annotated BayesC as a first-class manual topic without
overloading the existing single-trait annotated BayesC page.

## Current Problem

The current manual page at `docs/src/manual/annotated_bayesc.md` mixes two
different production paths:

- single-trait annotated BayesC
- dense 2-trait annotated BayesC

That makes the page harder to scan because the 2-trait path now has its own:

- 4-state joint prior over `00`, `10`, `01`, `11`
- 3-step tree annotation parameterization
- startup `Pi` semantics
- explicit `multi_trait_sampler=:auto|:I|:II` option
- support restrictions

## Design Decision

Keep `docs/src/manual/annotated_bayesc.md` as the single-trait-first guide and
add a separate page:

- `docs/src/manual/multitrait_annotated_bayesc.md`

The new page will be linked from:

- `docs/src/manual/annotated_bayesc.md`
- `docs/src/manual/workflow.md`
- `docs/make.jl`

## Scope Of The New Page

The new manual page should explain:

1. what multi-trait annotated BayesC is in JWAS today
2. the four joint marker states `00`, `10`, `01`, `11`
3. the 3-step tree prior:
   - `00` vs active
   - `11` vs singleton
   - `10` vs `01`
4. supported startup `Pi` behavior:
   - default `Pi=0.0` means all mass on `11` at startup
   - user-supplied joint `Pi` dictionary must sum to one
   - startup shared mass in `11` must be positive
5. supported runtime settings:
   - exactly 2 traits
   - dense only
   - `constraint=false`
   - `fast_blocks=false`
6. sampler selection:
   - `multi_trait_sampler=:auto`
   - `multi_trait_sampler=:I`
   - `multi_trait_sampler=:II`
7. output interpretation:
   - mean joint `pi_<geno>`
   - annotation coefficient step names
8. one worked example

## Changes To Existing Single-Trait Page

`docs/src/manual/annotated_bayesc.md` should be simplified so it stays the
entry point for annotated BayesC generally, but the 2-trait material should be
reduced to:

- a short note that JWAS also supports dense 2-trait annotated BayesC
- a direct link to the new multi-trait page

Single-trait dense, block, and streaming examples should remain on the original
page.

## Navigation

The docs navigation should expose the new page directly under the Manual
section, near the existing annotated BayesC and annotated BayesR pages.

Recommended order:

- Annotated BayesC
- Multi-Trait Annotated BayesC
- Annotated BayesR

## Documentation Style

The new page should be practical and production-focused:

- state current support clearly
- distinguish defaults from optional settings
- avoid long theoretical digressions
- keep notation aligned with the implementation and current manual language

## Verification

After the doc edits:

- build the docs with `julia --project=docs --startup-file=no docs/make.jl`
- confirm the new page is included in navigation and internal links resolve
