# Annotated BayesR Benchmark Implementation

## Goal

Add a production-path benchmark for dense annotated BayesR and compare it
against dense BayesR and dense annotated BayesC on a synthetic dataset with
informative annotations.

## Files Changed

Created:

- `benchmarks/annotated_bayesr_comparison.jl`
- `benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md`
- `docs/plans/2026-03-28-annotated-bayesr-benchmark-design.md`
- `docs/plans/2026-03-28-annotated-bayesr-benchmark-plan.md`
- `docs/plans/2026-03-28-annotated-bayesr-benchmark-implementation.md`

Modified:

- `test/unit/test_bayesr_parity.jl`

## What Was Implemented

### Benchmark script

`benchmarks/annotated_bayesr_comparison.jl` now:

- generates a synthetic genotype matrix and phenotype
- generates two marker-annotation columns
- builds truth metadata for:
  - true class
  - causal vs noncausal
  - enriched vs baseline annotation status
- runs production JWAS for:
  - `BayesR`
  - `Annotated_BayesR`
  - `Annotated_BayesC`
- writes:
  - `comparison_runs.csv`
  - `comparison_summary.csv`
  - `pip_group_summary.csv`
  - `annotation_coefficients.csv`
  - `truth_metadata.csv`

### Smoke test

`test/unit/test_bayesr_parity.jl` now includes a subprocess test that:

- runs the new benchmark script on a small synthetic setting
- checks that the expected output files are written
- checks the core schema of the summary outputs

## TDD Record

### Red

Added the subprocess smoke test first, then ran:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

It failed because:

- `benchmarks/annotated_bayesr_comparison.jl` did not exist
- expected CSV outputs were missing

### Green

Implemented the benchmark script, then ran a direct smoke command:

```bash
JWAS_ANNOT_BENCH_CHAIN_LENGTH=60 \
JWAS_ANNOT_BENCH_BURNIN=20 \
JWAS_ANNOT_BENCH_N_OBS=30 \
JWAS_ANNOT_BENCH_N_MARKERS=40 \
JWAS_ANNOT_BENCH_SEEDS=2026,2027 \
julia --project=. --startup-file=no benchmarks/annotated_bayesr_comparison.jl /tmp/annotated_bayesr_bench_smoke2
```

Then re-ran the driving test:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

and the benchmark smoke test passed.

## Full Benchmark Run

Full production benchmark command:

```bash
julia --project=. --startup-file=no benchmarks/annotated_bayesr_comparison.jl /tmp/annotated_bayesr_benchmark_20260328
```

Full output directory:

- `/tmp/annotated_bayesr_benchmark_20260328`

Settings used:

- `n_obs = 200`
- `n_markers = 1000`
- `chain_length = 10000`
- `burnin = 2000`
- seeds `2026, 2027, 2028, 2029, 2030`

## Main Findings

From the full benchmark:

- annotated BayesR learned the informative annotation direction at step 1
- annotated BayesR improved enriched-vs-baseline and causal-vs-noncausal PIP
  separation relative to ordinary BayesR
- annotated BayesR improved true-class ordering of PIP
- annotated BayesR did not improve `cor(y, EBV)` on this sparse scenario
- the likely limiting factor is weak identification in the step-2 and step-3
  conditional annotation models because classes `3` and `4` were rare

## Report

The detailed benchmark interpretation is saved in:

- `benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md`
