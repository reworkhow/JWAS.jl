# Annotated BayesR Less-Sparse Benchmark Design

## Goal

Add a second production benchmark scenario for annotated BayesR where the upper
BayesR classes are less sparse, then compare it against the first sparse
scenario.

## Why This Extension

The first benchmark showed:

- annotated BayesR learned the informative annotation direction
- annotated BayesR improved prioritization and class ordering
- annotated BayesR hurt `cor(y, EBV)` on a sparse scenario with only `27`
  markers in classes `3` and `4` combined

That result is informative, but it does not tell us whether the weak prediction
performance is a general method problem or a sparse-regime problem. The next
useful test is a scenario where the upper conditional annotation models have
more data.

## Scenario Design

Add a benchmark scenario selector to
`benchmarks/annotated_bayesr_comparison.jl`.

Scenarios:

- `sparse_upper_classes`
  - current default
- `less_sparse_upper_classes`
  - new scenario

For the new scenario:

- keep the same:
  - `n_obs`
  - `n_markers`
  - annotation structure
  - target heritability
- change the true mixture probabilities so classes `2`, `3`, and `4` are less
  rare, especially among enriched markers

Recommended class probabilities:

- baseline markers:
  - `[0.90, 0.05, 0.03, 0.02]`
- enriched markers:
  - `[0.65, 0.15, 0.12, 0.08]`

This keeps the same qualitative annotation story:

- `annotation_1` increases nonzero rate
- `annotation_1` also increases upper-class rate

but gives the step-2 and step-3 conditional models much more information.

## Outputs

The benchmark outputs should include the scenario explicitly:

- add `scenario` column to:
  - `comparison_runs.csv`
  - `comparison_summary.csv`
  - `pip_group_summary.csv`
  - `annotation_coefficients.csv`
  - `truth_metadata.csv`

That makes the report comparison straightforward and keeps reruns unambiguous.

## Testing

Update the subprocess smoke test so it runs the benchmark with:

- `JWAS_ANNOT_BENCH_SCENARIO=less_sparse_upper_classes`

and verifies:

- the script succeeds
- the output files exist
- the summary tables include the `scenario` column
- the written scenario value is `less_sparse_upper_classes`

## Report Update

Update the existing benchmark report rather than creating a separate report
file. The report should become a two-scenario comparison:

- Scenario 1: sparse upper classes
- Scenario 2: less sparse upper classes

The key comparison question is:

- does annotated BayesR become more competitive on `cor(y, EBV)` when the upper
  conditional annotation models are better identified?
