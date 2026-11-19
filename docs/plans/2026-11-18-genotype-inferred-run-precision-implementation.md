# Genotype-inferred run precision implementation

Base: PR #175 revision `fc287c26027f3ade2637cc234873cfdc2510460c`.
Design: [Genotype-inferred run precision](2026-11-18-genotype-inferred-run-precision-design.md).

## Changes

- Added internal `infer_run_double_precision(mme)`, inspecting only each genotype
  term's element type. It does not scan genotype values or allocate a converted
  matrix.
- Made the helper the lazy default for runMCMC's existing double_precision
  keyword. An omitted keyword resolves to false for uniform Float32 storage and
  true for uniform Float64 storage, before the shared-genotype default is chosen.
- Preserved false for models without genotype terms and Float32 streaming
  placeholders. Mixed or unsupported genotype element types fail before run
  setup with an actionable message.
- Preserved explicit true/false behavior. Explicit true can still request a
  guarded Float64 conversion; explicit false does not downcast existing Float64
  genotype storage. No new public keyword or precision sentinel was introduced.
- Preserved no-copy import defaults, shared views, the current sampler arithmetic,
  priors and sampling schedules. Updated the API documentation, manual example
  and changelog so precision need only be selected at genotype import.
- Extended the test fixture with an opt-in omitted-precision path. Its existing
  default still passes an explicit matching precision, preserving prior tests
  and reproduction controls.

## Verification

Tests were added before production edits. The genuine RED run reported 270
passes, 9 failures and 12 errors: omitted Float64 runs selected the wrong mode,
mixed terms did not fail before setup, and the resolver did not yet exist.
An initial test assumption about the untouched MCMCinfo sentinel was corrected
to compare its captured initial identity; production code did not change for it.

The final inferred-precision suite passed 573/573 checks; the previous conversion
guard and QC regression suite also passed 57/57, for 630 focused checks. Twelve
paired production comparisons cover single-trait and multi-trait samplers I/II,
Float32/Float64, and seeds 73/123. Inferred and explicitly matched settings give
exactly matching exported EBV, PEV, marker mean, SD and inclusion frequency.
The fixtures check expected counts, unique key sets and validated animal/marker
joins, and joined samples were inspected. They also verify parent/block aliases
and the resulting model-state element types. These are finite-chain equivalence
checks, not convergence or prediction-accuracy studies.

The full default suite passed **2,397/2,397 checks** under Julia 1.11.7 on Apple
Silicon. Optional integration tests remained disabled by the default test
configuration. The documentation build passed against this worktree's source;
deployment was skipped. Independent scope and quality reviews found no blocking
issues. An additional independent production smoke audit passed 15/15 checks
for omitted precision with Float64 GBLUP, RR-BLUP and two BayesC genotype terms.

Validation commands:

```sh
GKSwstype=100 JULIA_PKG_OFFLINE=true julia --project=. --startup-file=no --compiled-modules=existing test/runtests.jl
GKSwstype=100 JULIA_PKG_OFFLINE=true JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=docs --startup-file=no --compiled-modules=existing docs/make.jl
```

The ignored docs environment points at this worktree, with a writable temporary
first Julia depot for runtime logs. Dependencies and package version are unchanged.

## Scope and historical measurements

The earlier RSS report under
`benchmarks/reports/2026-11-18-bayesc-no-copy-revision/` remains unchanged and
identifies its original measured source commit and fingerprint. Its 36-case RSS
grid was not rerun or relabelled as measurements of this inference follow-up.
No target-scale execution or new performance/convergence claim is made here.

Dense CSV/DataFrame intermediate allocations, actual SNP-filtering copies and
any redesign of storage versus computation precision remain separate work.
