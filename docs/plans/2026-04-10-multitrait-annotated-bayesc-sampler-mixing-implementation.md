# Multi-Trait Annotated BayesC Sampler Mixing Benchmark Implementation

## Summary

Added a focused follow-up benchmark for the sampler-I versus sampler-II
question in multi-trait annotated BayesC.

This follow-up keeps the earlier full benchmark as context, but narrows the
longer-chain rerun to the three multi-trait production methods:

- conventional multi-trait BayesC
- multi-trait annotated BayesC with `multi_trait_sampler = :I`
- multi-trait annotated BayesC with `multi_trait_sampler = :II`

It also adds an exact tiny-case regression so the sampler comparison is framed
correctly as a finite-chain mixing question, not a different-posterior
question.

## Files Changed

- Modified `benchmarks/simulated_annotations_multitrait_comparison.jl`
- Modified `test/unit/test_misc_coverage.jl`
- Modified `test/unit/test_multitrait_mcmc.jl`
- Added planning records:
  - `docs/plans/2026-04-10-multitrait-annotated-bayesc-sampler-mixing-design.md`
  - `docs/plans/2026-04-10-multitrait-annotated-bayesc-sampler-mixing-plan.md`
- Added report:
  - `benchmarks/reports/2026-04-10-multitrait-annotated-bayesc-sampler-mixing-report.md`

## What Changed

### Exact tiny-case regression

Extended `test/unit/test_multitrait_mcmc.jl` with a one-marker, two-trait
annotated BayesC check that:

- computes the exact four-state posterior over `00`, `10`, `01`, and `11`
- runs the internal sampler-I kernel
- runs the internal sampler-II kernel
- verifies both empirical state-frequency vectors stay within Monte Carlo
  tolerance of the same exact target

This regression is the correctness backstop for the mixing benchmark.

### Focused multi-trait rerun mode

Extended `benchmarks/simulated_annotations_multitrait_comparison.jl` with an
optional focused mode controlled by:

- `JWAS_SIMULATED_MT_FOCUSED_MULTITRAIT=true`

When enabled, the harness only runs:

- `MT_BayesC`
- `MT_Annotated_BayesC_I`
- `MT_Annotated_BayesC_II`

The default full benchmark matrix remains unchanged when the flag is absent.

### Shared-state mixing summaries

Added a compact multi-trait mixing summary built from the existing keyed
benchmark outputs.

The follow-up now reports, across seeds:

- shared precision, recall, and F1 from `P11` top-`k`
- mean `P11` on true shared loci
- mean `P11` on nonshared loci
- runtime
- variability across seeds for the shared metrics

The benchmark continues to reconstruct `P11` from saved MCMC samples and keeps
all marker/truth comparisons keyed by explicit marker identifier.

### Contract coverage

Extended `test/unit/test_misc_coverage.jl` so the benchmark contract now
checks:

- focused multi-trait case selection
- the sampler metadata for annotated `:I` and `:II` cases
- the shared posterior row contract
- the mixing-summary schema
- the focused-mode wiring through `main()`

## Benchmark Runs

### Smoke rerun

```bash
JWAS_SIMULATED_MT_FOCUSED_MULTITRAIT=true \
JWAS_SIMULATED_MT_SEEDS=11 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=200 \
JWAS_SIMULATED_MT_BURNIN=50 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=50 \
JWAS_SIMULATED_MT_WARMUP=false \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_sampler_mixing_smoke
```

### Production rerun

```bash
JWAS_SIMULATED_MT_FOCUSED_MULTITRAIT=true \
JWAS_SIMULATED_MT_SEEDS=101,202,303,404 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=4000 \
JWAS_SIMULATED_MT_BURNIN=1000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=40 \
JWAS_SIMULATED_MT_WARMUP=false \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_sampler_mixing_20260410
```

## Verification

Ran:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'`
- the smoke rerun above
- the production rerun above

## Result

The repository now has:

- an exact-regression backstop showing sampler I and sampler II target the
  same tiny-case posterior
- a focused production rerun that measures longer-chain finite-sample
  differences in shared-state recovery and runtime
- a follow-up report that interprets sampler-I versus sampler-II as a mixing
  and efficiency comparison rather than a different-model comparison
