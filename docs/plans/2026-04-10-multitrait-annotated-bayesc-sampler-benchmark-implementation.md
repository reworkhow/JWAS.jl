# Multi-Trait Annotated BayesC Sampler Benchmark Implementation

## Summary

Extended the production benchmark harness for the packaged
`simulated_annotations` 2-trait fixture so it now compares:

- conventional multi-trait BayesC
- multi-trait annotated BayesC with Gibbs sampler I
- multi-trait annotated BayesC with Gibbs sampler II
- single-trait BayesC
- single-trait annotated BayesC
- single-trait BayesR
- single-trait annotated BayesR

The benchmark now reports true `P11`-based pleiotropic recovery for the
multi-trait methods and overlap-based shared detection for the paired
single-trait baselines.

## Files Changed

- Modified `benchmarks/simulated_annotations_multitrait_comparison.jl`
- Modified `test/unit/test_misc_coverage.jl`
- Added planning records:
  - `docs/plans/2026-04-10-multitrait-annotated-bayesc-sampler-benchmark-design.md`
  - `docs/plans/2026-04-10-multitrait-annotated-bayesc-sampler-benchmark-plan.md`
- Added benchmark report:
  - `benchmarks/reports/2026-04-10-multitrait-annotated-bayesc-sampler-benchmark-report.md`

## What Changed

### Benchmark method matrix

Expanded the method matrix in
`benchmarks/simulated_annotations_multitrait_comparison.jl` to include:

- `MT_BayesC`
- `MT_Annotated_BayesC_I`
- `MT_Annotated_BayesC_II`
- `BayesC_y1`, `BayesC_y2`
- `Annotated_BayesC_y1`, `Annotated_BayesC_y2`
- `BayesR_y1`, `BayesR_y2`
- `Annotated_BayesR_y1`, `Annotated_BayesR_y2`

The multi-trait annotated cases now pass `multi_trait_sampler=:I` and `:II`
through the production `get_genotypes` path.

### Shared-recovery summaries

Reworked the shared-recovery logic so the benchmark now distinguishes the
multi-trait and single-trait cases explicitly.

For multi-trait methods:

- reconstruct `P11 = Pr(delta = (1,1) | data)` from the saved marker-effect
  MCMC traces
- rank markers by `P11`
- declare the top `k_shared_true` markers as shared
- summarize shared precision, recall, F1, and `P11` separation

For single-trait method families:

- pair the trait-1 and trait-2 joined marker tables by normalized `marker_id`
- detect top `k_y1_true` and top `k_y2_true` markers separately
- declare shared markers by overlap of those detected sets
- summarize shared precision, recall, and F1 for the overlap proxy

### Contract coverage

Strengthened `test/unit/test_misc_coverage.jl` so it now checks:

- the actual method-matrix metadata returned by `method_cases()`
- that only the annotated multi-trait BayesC cases carry non-`auto`
  `multi_trait_sampler` settings
- that `shared_posterior_row` uses the `P11_topk` contract
- that `single_trait_shared_overlap_row` uses the overlap-based
  single-trait contract
- that `main()` wires the benchmark through
  `summarize_multitrait_shared_posteriors` and
  `summarize_single_trait_shared_overlaps`

## Benchmark Runs

### Smoke run

```bash
JWAS_SIMULATED_MT_SEEDS=11 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=200 \
JWAS_SIMULATED_MT_BURNIN=50 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=50 \
JWAS_SIMULATED_MT_WARMUP=false \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_sampler_benchmark_smoke
```

### Production run

```bash
JWAS_SIMULATED_MT_SEEDS=101,202 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=2000 \
JWAS_SIMULATED_MT_BURNIN=500 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=20 \
JWAS_SIMULATED_MT_WARMUP=false \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_sampler_benchmark_20260410
```

## Verification

Run:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'`
- `julia --project=. --startup-file=no test/runtests.jl`
- the smoke benchmark command above
- the production benchmark command above

## Result

The benchmark harness now supports the intended sampler-I versus sampler-II
comparison for multi-trait annotated BayesC and produces a reportable
production comparison against the requested single-trait BayesC/BayesR
baselines on the packaged simulated annotation fixture.
