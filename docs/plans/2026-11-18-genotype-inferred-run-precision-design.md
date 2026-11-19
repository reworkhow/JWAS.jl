# Genotype-inferred run precision design

Base: PR #175 revision `fc287c26027f3ade2637cc234873cfdc2510460c`.

## Objective

Let runMCMC inherit the genotype matrix's precision when double_precision is
omitted. Users should select genotype precision once at import, avoiding a second
setting that can cause a whole-matrix promotion. Preserve no-copy adoption,
shared views, existing explicitly requested conversion and its memory guard.

## Interface and behavior

- Use a lazy keyword default resolver, returning a Bool before any later keyword
  defaults inspect double_precision. No sentinel is stored in MCMCinfo and no new
  public precision keyword is introduced.
- All genotype terms Float32: use double_precision=false. All Float64: use true.
  Inspect element types only; precision inference must not traverse/copy data.
- No genotype terms: preserve the conventional-model false default. Streaming
  genotypes remain Float32 and must not be promoted by inference.
- Mixed Float32/Float64 terms: fail with an actionable message before run setup
  and require an explicit choice. Never silently promote a large genotype term
  solely to reconcile an inferred global precision.
- Explicit true/false keep existing behavior, including guarded Float64 promotion
  and strict shared-storage validation. Priors, random draws and sampler arithmetic
  are not otherwise changed.
- Automatic shared-genotype eligibility sees the resolved Bool, so an imported
  Float64 BayesC matrix can keep its original parent and shared rows without
  repeating double_precision=true in runMCMC.

## Test-first implementation and verification

Add production regression tests before changing code: omitted precision with
Float32/Float64 input, parent identity and state types, conventional fallback,
mixed-term early error, explicit overrides, and streaming Float32 resolution.
Compare automatic and explicitly matched runs with actual production sampling,
explicit animal/marker keys, expected counts and unique matching key sets.
Use multiple seeds where stochastic comparisons are made.

Run focused tests and independently review the resolver and keyword-default
ordering. Then run the full default test suite and documentation build against
this worktree. Document inferred versus explicit precision and record the exact
validation in a matching implementation note before updating the existing PR.

## Scope limits

This adds a precision-selection policy, not a new numerical accumulation scheme.
Existing mixed arithmetic remains. CSV loading, actual marker filtering and the
broader storage/computation precision API remain separate. Earlier benchmark
reports retain their original source fingerprints; do not relabel measurements
from a prior revision as measurements of this one.
