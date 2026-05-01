# Release Polish Documentation Implementation

**Date:** 2026-05-01

**Status:** completed

## Summary

This docs-only release polish pass updated the user-facing release surface for
recent annotated, streaming, and fast-block work:

- `CHANGELOG.md` now has a 2026-05-01 release note section covering
  independent-block fast-block documentation, dense 2-trait annotated BayesC,
  simulated annotation datasets, streaming genotype workflows, and support
  boundaries.
- `README.md` now describes JWAS as a Bayesian mixed-model package for genomic
  prediction and GWAS with single-trait and multi-trait analyses, dense,
  streaming, fast-block genotype workflows, and annotated BayesC/BayesR support.
  It now points current users to the Documenter manual first while retaining the
  wiki link for legacy and community examples.
- `docs/src/index.md` now makes annotated BayesC/BayesR, dense 2-trait
  annotated BayesC, streaming genotype storage, exact fast blocks, and
  approximate independent blocks easier to discover from the landing page.
- Annotated method pages were aligned across
  `docs/src/manual/annotated_bayesc.md`,
  `docs/src/manual/annotated_bayesr.md`, and
  `docs/src/manual/multitrait_annotated_bayesc.md`, with support boundaries,
  sampler choices, streaming constraints, and dense 2-trait guidance kept near
  the relevant examples.
- `docs/src/manual/public.md` now includes the current public entry point for
  `outputMCMCsamples`.
- `docs/src/manual/block_bayesc.md` now more directly distinguishes exact
  sequential fast-block sweeps from the opt-in approximate
  `independent_blocks=true` schedule, including explicit block-start semantics
  and the `OPENBLAS_NUM_THREADS=1` guidance for threaded block runs.

No Julia source, benchmark scripts, generated assets, README/changelog/manual
pages outside the release-polish scope, or unrelated files were changed during
Task 5.

## Wording Decisions

- Kept `independent_blocks=false` as the documented default exact sequential
  fast-block sweep.
- Described `independent_blocks=true` as opt-in and approximate unless off-block
  weighted genotype crossproducts are effectively zero.
- Kept the README examples path docs-first, with the wiki framed as legacy and
  community material rather than the primary source of current examples.
- Avoided adding single-step or old `axpy!` debugging material because it was
  out of scope for this docs-only release polish pass.
- Avoided new benchmark claims; existing performance-oriented language remains
  tied to the documented workflow descriptions.

## Verification

Ran:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Result: exit code 0. Documenter emitted the acceptable deployment environment
auto-detection warning and skipped deployment.

Ran:

```bash
rg -n "singe-step|old wiki|TODO|FIXME|independent block|independent_blocks|storage=:stream" README.md CHANGELOG.md docs/src
```

Result: expected release-polish hits remain for `independent_blocks`,
`independent block`, and `storage=:stream`. No `singe-step`, `old wiki`,
`TODO`, or `FIXME` hits were present in the searched files.

## Benchmarks And Artifacts

Production benchmarks were not rerun for this docs polish task. Unrelated
untracked benchmark and debug artifacts were not touched.
