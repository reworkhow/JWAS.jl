# Multi-Trait Annotated BayesC Manual Implementation

## Summary

Split the dense 2-trait annotated BayesC documentation into its own manual
page and rewired the existing annotated BayesC documentation to keep the
single-trait page as the main entry point.

## Files Changed

- Created `docs/src/manual/multitrait_annotated_bayesc.md`
- Updated `docs/src/manual/annotated_bayesc.md`
- Updated `docs/src/manual/workflow.md`
- Updated `docs/make.jl`
- Added planning records:
  - `docs/plans/2026-04-09-multitrait-annotated-bayesc-manual-design.md`
  - `docs/plans/2026-04-09-multitrait-annotated-bayesc-manual-plan.md`

## What Changed

### New dedicated manual page

Added a dedicated page for dense 2-trait annotated BayesC that documents:

- the 4-state prior over `00`, `10`, `01`, `11`
- the 3-step tree annotation parameterization
- startup `Pi` defaults and validation requirements
- explicit `multi_trait_sampler=:auto|:I|:II`
- current support restrictions
- a worked dense 2-trait example
- output interpretation for the joint `pi_<geno>` summary and annotation steps

### Existing annotated BayesC page

Reframed the original `annotated_bayesc.md` page so it stays focused on the
single-trait dense, block, and streaming paths, with a short dedicated section
that points readers to the new multi-trait page.

### Navigation and workflow

Added the new page to the Documenter manual navigation and linked it from the
general workflow page where annotated BayesC methods are introduced.

## Verification

Built the docs successfully with:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

The build completed successfully. The only output note was the usual
Documenter warning that deployment was skipped because no deployment
environment was detected.
