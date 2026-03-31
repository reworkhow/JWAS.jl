# Annotated BayesR Less-Sparse Benchmark Implementation

## Goal

Extend the production annotated BayesR benchmark with a second synthetic
scenario where upper BayesR classes are less sparse, then update the benchmark
report to compare both scenarios clearly.

## Changes

### Benchmark script

Updated:

- `benchmarks/annotated_bayesr_comparison.jl`

Changes made:

- added scenario selection through `JWAS_ANNOT_BENCH_SCENARIO`
- kept `sparse_upper_classes` as the default
- added `less_sparse_upper_classes` with:
  - baseline class probabilities `[0.90, 0.05, 0.03, 0.02]`
  - enriched class probabilities `[0.65, 0.15, 0.12, 0.08]`
- carried `scenario` through:
  - `truth_metadata.csv`
  - `comparison_runs.csv`
  - `comparison_summary.csv`
  - `pip_group_summary.csv`
  - `annotation_coefficients.csv`

### Tests

Updated:

- `test/unit/test_bayesr_parity.jl`

Changes made:

- extended the annotated BayesR benchmark smoke test to run with
  `JWAS_ANNOT_BENCH_SCENARIO=less_sparse_upper_classes`
- added assertions that:
  - `comparison_summary.csv` contains a `scenario` column
  - `truth_metadata.csv` contains a `scenario` column
  - both files carry the expected scenario label

### Report

Updated:

- `benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md`

Changes made:

- converted the original single-scenario report into a two-scenario comparison
- added the new less-sparse truth summary
- added scenario-specific method tables
- updated the interpretation to separate:
  - mechanical correctness
  - prioritization behavior
  - prediction behavior
  - upper-step identification issues

## Benchmark Runs

Existing sparse-scenario run used for comparison:

- `/tmp/annotated_bayesr_benchmark_20260328`

New less-sparse run:

- `/tmp/annotated_bayesr_benchmark_20260328_less_sparse`

Command used for the new run:

```bash
JWAS_ANNOT_BENCH_SCENARIO=less_sparse_upper_classes \
  julia --project=. --startup-file=no \
  benchmarks/annotated_bayesr_comparison.jl \
  /tmp/annotated_bayesr_benchmark_20260328_less_sparse
```

## Verification

Targeted benchmark smoke test:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Planned full-branch verification after the report update:

```bash
julia --project=. --startup-file=no test/runtests.jl
julia --project=docs --startup-file=no docs/make.jl
```

## Conclusion

The benchmark extension is small in code surface but materially improves the
scientific interpretation:

- the sparse scenario alone made annotated BayesR look more fragile than it is
- the less-sparse scenario shows the same prioritization gains with a smaller
  prediction penalty
- the persistent weakness is still in the upper conditional annotation models,
  not in the basic production implementation
