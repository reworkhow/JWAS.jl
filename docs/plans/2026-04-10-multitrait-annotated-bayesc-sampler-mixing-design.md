# Multi-Trait Annotated BayesC Sampler Mixing Benchmark Design

## Goal

Follow up the initial sampler benchmark with a more careful production-path
comparison focused on the remaining method question:

- do sampler I and sampler II differ because they target different posteriors,
  or because they mix differently at the finite chain lengths used in the first
  benchmark?

The new exact one-marker regression test already establishes that both kernels
target the same posterior on a tiny case. This benchmark extension should now
focus on practical mixing and efficiency on the packaged
`simulated_annotations` two-trait fixture.

## Scope

Keep the production benchmark on the packaged fixture under
`src/4.Datasets/data/simulated_annotations`, but narrow the longer-chain rerun
to the multi-trait methods that matter for the sampler question:

- conventional multi-trait BayesC
- multi-trait annotated BayesC with `multi_trait_sampler=:I`
- multi-trait annotated BayesC with `multi_trait_sampler=:II`

Single-trait BayesC and BayesR baselines should remain in the existing report
for context, but they do not need the same longer-chain rerun unless a later
question requires it.

## Benchmark Question

The new report should answer:

1. whether sampler I and sampler II still show materially different
   pleiotropic recovery after using longer chains
2. whether the difference is accompanied by lower between-seed variability for
   sampler II
3. what runtime cost sampler II pays for any shared-state mixing advantage

This is a mixing-efficiency comparison, not a correctness comparison.

## Metrics

### Primary

Primary emphasis remains pleiotropic recovery, using the true joint posterior
shared-state score:

- `P11 = Pr(delta = (1,1) | data)`

For each seed and method:

- rank markers by `P11`
- declare the top `k_shared_true` markers as shared
- report shared precision, recall, and F1

### Diagnostics

Add diagnostics that help interpret finite-chain differences:

- mean `P11` on true shared loci
- mean `P11` on true nonshared loci
- per-seed shared precision, recall, and F1
- summarized variability across seeds for those shared metrics
- runtime

The per-seed summaries are important. If sampler II still looks stronger after
longer chains, and does so with lower variability, that supports a mixing
advantage rather than a benchmark artifact.

### Secondary

Keep:

- any-active recovery
- per-trait top-k recall
- genomic prediction accuracy `cor(y, EBV)`

These remain secondary because the current question is specifically about the
multi-trait sampler behavior on pleiotropic recovery.

## Execution

Use the same production benchmark harness in
`benchmarks/simulated_annotations_multitrait_comparison.jl`.

Recommended execution:

1. keep the current benchmark outputs and report as the initial comparison
2. add an optional longer-chain focused rerun mode for the three multi-trait
   methods
3. write a follow-up report that compares:
   - conventional multi-trait BayesC
   - annotated sampler I
   - annotated sampler II
4. interpret results in light of the new exact posterior regression test

## Comparison Rules

Preserve the repository result-comparison rules:

- join markers to truth by explicit marker key
- normalize marker identifiers before joining
- assert key-set equality and join row counts
- save comparison-ready summaries with explicit labels

Do not infer anything from row order.

## Expected Outcome

The useful outputs of this follow-up are:

- a benchmark run with longer chains on the three multi-trait methods
- per-seed and summarized shared-state diagnostics
- a report that clearly separates:
  - same target posterior
  - different finite-chain mixing behavior

That will make the sampler-I versus sampler-II tradeoff defensible without
claiming the kernels represent different models.
