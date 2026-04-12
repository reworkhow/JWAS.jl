# Simulated Annotations Cross-Validation Implementation

## Scope

Implemented production K-fold cross-validation benchmarking on the packaged
`simulated_annotations` fixture inside
`benchmarks/simulated_annotations_multitrait_comparison.jl`.

The implementation adds:

- deterministic fold assignment by `ID` and seed
- phenotype masking for held-out individuals
- held-out EBV scoring for multi-trait and single-trait methods
- CV-specific method and family summaries
- support for the full current multi-trait model set, including the
  empty-annotated controls

## Files Changed

- modified:
  - `benchmarks/simulated_annotations_multitrait_comparison.jl`
  - `test/unit/test_misc_coverage.jl`
- added:
  - `benchmarks/reports/2026-04-11-simulated-annotations-cv-report.md`
  - `docs/plans/2026-04-11-simulated-annotations-cv-design.md`
  - `docs/plans/2026-04-11-simulated-annotations-cv-plan.md`

## Commands Run

Targeted contract test:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Smoke CV benchmark:

```bash
JWAS_SIMULATED_MT_BENCHMARK_MODE=cv \
JWAS_SIMULATED_MT_CV_FOLDS=2 \
JWAS_SIMULATED_MT_SEEDS=11 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=150 \
JWAS_SIMULATED_MT_BURNIN=50 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=50 \
JWAS_SIMULATED_MT_WARMUP=false \
julia --project=. --startup-file=no benchmarks/simulated_annotations_multitrait_comparison.jl /tmp/jwas_simulated_annotations_cv_smoke2
```

Production CV benchmark:

```bash
JWAS_SIMULATED_MT_BENCHMARK_MODE=cv \
JWAS_SIMULATED_MT_CV_FOLDS=5 \
JWAS_SIMULATED_MT_SEEDS=101,202 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=1500 \
JWAS_SIMULATED_MT_BURNIN=500 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=50 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no benchmarks/simulated_annotations_multitrait_comparison.jl /tmp/jwas_simulated_annotations_cv_full_20260411
```

Full test suite:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

## Output Files

Main production output directory:

- `/tmp/jwas_simulated_annotations_cv_full_20260411`

Key generated files:

- `cv_fold_assignments.csv`
- `cv_per_fold_summary.csv`
- `cv_method_summary.csv`
- `cv_family_summary.csv`

## Main Results

See the benchmark report for the full interpretation:

- `benchmarks/reports/2026-04-11-simulated-annotations-cv-report.md`

High-level results:

- `MT_Annotated_BayesC_II` had the best multi-trait family mean held-out
  correlation
- `BayesR_single` had the best single-trait family mean held-out correlation
- `MT_BayesC` and `MT_BayesC_I` were identical on held-out prediction in this
  benchmark, as expected from `:auto` dispatch

## Notes

- Multi-trait CV hides both traits in the held-out fold. This keeps the
  comparison aligned with the requested single-trait baselines.
- Empty-annotated multi-trait BayesC remains a diagnostic control. It is
  included in the CV run because the requested multi-trait comparison covered
  all current multi-trait models.
