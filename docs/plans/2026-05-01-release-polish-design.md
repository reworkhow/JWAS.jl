# Release Polish Documentation Design

**Date:** 2026-05-01

**Status:** approved

## Goal

Polish JWAS user-facing documentation around the recent production updates so
users can understand the current feature set without reading internal plan
files or old issue threads.

This is a documentation polish pass. It should not change sampler behavior,
benchmark scripts, or production code.

## Scope

The polish pass should focus on recent JWAS work:

- annotated BayesC and annotated BayesR
- dense 2-trait annotated BayesC
- streaming genotype workflows
- exact fast-block sampling
- approximate independent-block fast-block sampling
- release-note clarity for the current public package state

Single-step debugging and the older fixed `axpy!` issue are out of scope for
this polish pass. They are separate maintenance work and should not distract
from the recent genomic sampler and documentation updates.

## User-Facing Story

The revised docs should present JWAS as having three production-oriented
genomic-data paths:

1. Dense genotype workflows for standard analyses and annotated methods.
2. Streaming genotype workflows for memory-sensitive single-trait BayesC use.
3. Fast-block workflows for larger dense genotype analyses.

Within fast-block workflows, the docs should distinguish:

- `fast_blocks` with `independent_blocks=false`: exact sequential blocked
  sweep.
- `independent_blocks=true`: approximate independent-block sweep that can use
  block-level parallelism when the user accepts the off-block approximation.

This distinction should be repeated wherever users are likely to encounter the
option.

## Files To Consider

Expected edits:

- `CHANGELOG.md`: summarize recent user-facing changes and docs additions.
- `README.md`: lightly refresh the project overview and direct users to the
  current manual rather than old wiki-first workflows.
- `docs/src/index.md`: improve discoverability of annotated, streaming, and
  fast-block topics.
- `docs/src/manual/public.md`: ensure current public functions and options are
  easy to find.
- `docs/src/manual/annotated_bayesc.md`: tighten single-trait guidance and link
  clearly to multi-trait annotated BayesC.
- `docs/src/manual/annotated_bayesr.md`: align wording with the current
  annotated-method story.
- `docs/src/manual/multitrait_annotated_bayesc.md`: ensure support boundaries
  and sampler choices are clear.
- `docs/src/manual/block_bayesc.md`: clarify exact blocked sweeps versus
  independent-block approximation.
- Streaming manual pages: make only small consistency edits if needed.

The implementation should stay concise. The goal is a clearer release surface,
not a full manual rewrite.

## Caveats To Document

The docs should call out these constraints where relevant:

- `independent_blocks=true` requires `fast_blocks != false`.
- Independent blocks are approximate unless off-block weighted genotype
  crossproducts are effectively zero.
- User-provided block starts use full-sweep chain semantics.
- Streaming genotype storage is memory-oriented and has precision/storage
  constraints.
- Annotated methods require valid annotation matrices and supported method
  combinations.
- Dense 2-trait annotated BayesC has explicit support boundaries that should
  remain visible near examples.

## Verification

Verification should be documentation-focused:

- run `julia --project=docs --startup-file=no docs/make.jl`
- inspect changed links and navigation entries
- run a quick text search for stale or contradictory wording in changed pages

Do not run production benchmarks as part of this polish pass.

## Non-Goals

- no sampler changes
- no benchmark reruns
- no new public API
- no cleanup of unrelated benchmark/debug artifacts
- no single-step issue work in this polish branch
