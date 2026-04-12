# Multi-Trait BayesC Plain vs Annotated Sampler Benchmark Design

## Goal

Follow up the current sampler discussion with a production-path comparison that
answers a narrower question:

- does the sampler-I versus sampler-II finite-chain gap already appear in
  plain multi-trait BayesC, or does it mainly emerge once the prior becomes
  marker-specific under annotation?

The existing one-marker exact regression already supports a shared target
posterior for the annotated samplers. The next step should move to a more
complicated production dataset without introducing a new simulation design.

## Scope

Use the packaged two-trait `simulated_annotations` fixture under
`src/4.Datasets/data/simulated_annotations` as the production dataset.

Compare four explicit multi-trait production runs:

- `MT_BayesC_I`
- `MT_BayesC_II`
- `MT_Annotated_BayesC_I`
- `MT_Annotated_BayesC_II`

This keeps the dataset fixed and varies only:

- plain vs annotated prior structure
- sampler I vs sampler II

The older `MT_BayesC` auto-dispatch row can remain in the broader benchmark
for backward compatibility, but it should not be the primary row in this
focused follow-up because the question is now about explicit sampler choice.

## Benchmark Question

The follow-up should answer:

1. whether plain multi-trait BayesC also shows a practical sampler-I versus
   sampler-II gap on the packaged two-trait fixture
2. whether the annotation layer amplifies that gap
3. whether sampler-II runtime cost is similar in the plain and annotated
   settings

This is still a finite-chain mixing and efficiency comparison, not a claim that
the samplers target different posteriors.

## Metrics

### Primary

Keep pleiotropic recovery first.

For all four multi-trait methods, use the true joint posterior shared-state
score:

- `P11 = Pr(delta = (1,1) | data)`

For each seed and method:

- reconstruct `P11` from saved MCMC samples
- rank markers by `P11`
- declare the top `k_shared_true` markers as shared
- report shared precision, recall, and F1

### Secondary

Keep the same secondary metrics as the current multi-trait benchmark:

- any-active recovery
- per-trait top-k recall
- genomic prediction accuracy `cor(y, EBV)`
- runtime

### Diagnostics

Preserve:

- mean `P11` on true shared loci
- mean `P11` on nonshared loci
- per-seed shared metrics
- variability across seeds

Those are needed to separate “same target, different mixing” from simple
benchmark noise.

## Test Backstop

Add a plain multi-trait BayesC tiny-case regression alongside the existing
annotated one.

That regression should:

- reuse the exact-posterior helper
- switch the prior provider from `MarkerSpecificPiPrior` to `GlobalPiPrior`
- verify sampler I and sampler II both match the same exact posterior within
  Monte Carlo tolerance

This gives symmetric low-level coverage for:

- plain multi-trait BayesC
- annotated multi-trait BayesC

## Execution

Reuse the existing production benchmark harness in
`benchmarks/simulated_annotations_multitrait_comparison.jl`.

Recommended execution:

1. add explicit plain multi-trait BayesC sampler-I and sampler-II method cases
2. add a focused selector for the four-method sampler comparison
3. run a smoke benchmark to verify wiring
4. run a longer-chain production rerun on those four cases
5. write a follow-up report comparing plain and annotated sampler behavior

## Comparison Rules

Preserve the repository comparison rules:

- join markers by explicit marker key
- normalize marker IDs before joining
- assert join counts and key-set equality
- do not compare by row order

Do not weaken any of the existing benchmark alignment checks.

## Expected Outcome

The useful output of this follow-up is not just another sampler table. It is a
clear separation between:

- a **general** sampler-I versus sampler-II finite-chain effect that already
  exists in plain multi-trait BayesC
- an **annotation-specific** amplification, if the gap becomes larger once the
  prior is marker-specific

That will tell us whether the practical sampler discussion is mainly about the
multi-trait Gibbs kernels themselves or about their interaction with the
annotated prior.
