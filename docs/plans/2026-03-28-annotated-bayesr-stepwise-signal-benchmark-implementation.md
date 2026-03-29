# Annotated BayesR Stepwise-Signal Benchmark Implementation

Date: 2026-03-28

## Goal

Add a third production benchmark scenario for annotated BayesR where the
annotation truth directly drives all three sequential conditional probabilities:

- `P(delta > 1)`
- `P(delta > 2 | delta > 1)`
- `P(delta > 3 | delta > 2)`

The purpose was to test whether the step-2 and step-3 annotation models become
identifiable when the synthetic truth is expressed in the same sequential form
as the sampler.

## Changes

### 1. Benchmark script

Updated:

- `benchmarks/annotated_bayesr_comparison.jl`

Changes:

- added `joint_probs_from_conditionals(p1, p2, p3)`
- added scenario `stepwise_annotation_signal`
- defined the new truth through sequential conditional probabilities:
  - baseline:
    - `p1 = 0.10`
    - `p2 = 0.20`
    - `p3 = 0.20`
  - enriched:
    - `p1 = 0.30`
    - `p2 = 0.60`
    - `p3 = 0.60`
- converted those conditionals into 4-class joint probabilities inside the
  benchmark data generator

### 2. Smoke test

Updated:

- `test/unit/test_bayesr_parity.jl`

Changes:

- changed the annotated BayesR production benchmark smoke test to run the new
  scenario
- verified the new scenario name is propagated into:
  - `comparison_summary.csv`
  - `truth_metadata.csv`

### 3. Report

Updated:

- `benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md`

Changes:

- expanded the report from two scenarios to three
- added the new scenario truth definition and realized counts
- added the new production benchmark results
- updated the interpretation and bottom-line conclusions

## Benchmark Result

Output directory:

- `/tmp/annotated_bayesr_benchmark_20260328_stepwise_signal`

Main result:

- annotated BayesR moved the informative annotation in the correct direction
  for all three steps:
  - step 1 `Annotation_1 = +1.3349`
  - step 2 `Annotation_1 = +0.2365`
  - step 3 `Annotation_1 = +0.2518`
- annotated BayesR improved prioritization over ordinary BayesR:
  - top-causal recall `0.2986` vs `0.2222`
  - enriched minus baseline PIP gap `0.0963` vs `0.0067`
- but `cor(y, EBV)` remained below ordinary BayesR:
  - `0.7409` vs `0.9164`

Interpretation:

- the deeper sequential annotation models are responsive when the truth
  explicitly drives them
- the remaining concern is not basic mechanics
- the remaining concern is the tradeoff between improved prioritization and
  weaker prediction, together with nuisance-annotation movement

## Verification

Executed:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'`
- `JWAS_ANNOT_BENCH_SCENARIO=stepwise_annotation_signal julia --project=. --startup-file=no benchmarks/annotated_bayesr_comparison.jl /tmp/annotated_bayesr_benchmark_20260328_stepwise_signal`

Planned as final branch verification after these benchmark/report edits:

- `julia --project=. --startup-file=no test/runtests.jl`
- `julia --project=docs --startup-file=no docs/make.jl`
