# Simulated Annotations Multi-Trait Benchmark Implementation

## Summary

Extended the packaged `simulated_annotations` fixture with native 2-trait files,
added a production-path benchmark harness for dense multi-trait annotated
BayesC and the requested baselines, and recorded the benchmark outputs in a
dedicated report.

## Code and Data Changes

### Packaged fixture

Updated
[generate_dataset.R](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/4.Datasets/data/simulated_annotations/generate_dataset.R)
to generate additional 2-trait files:

- `phenotypes_mt.csv`
- `annotations_mt.csv`
- `truth_mt.csv`

Updated
[README.md](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/4.Datasets/data/simulated_annotations/README.md)
to document:

- the 4-state 2-trait truth
- the effect simulation
- the 3-step annotation signal design
- the new saved files

### Tests

Updated
[test_misc_coverage.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/test/unit/test_misc_coverage.jl)
to verify the new fixture files resolve through `JWAS.Datasets` and expose the
expected columns.

### Benchmark harness

Added
[simulated_annotations_multitrait_comparison.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/benchmarks/simulated_annotations_multitrait_comparison.jl),
which:

- resolves the packaged 2-trait fixture through `JWAS.Datasets`
- runs dense production-path:
  - multi-trait BayesC
  - multi-trait annotated BayesC
  - single-trait BayesC / annotated BayesC on each trait
  - single-trait BayesR on each trait
- enforces marker and ID joins by explicit normalized keys
- writes joined marker tables, per-run summaries, aggregate method summaries,
  and annotation-coefficient summaries
- performs a small unrecorded warmup per case so runtime rows are not dominated
  by first-call compilation

## Verification

Run:

- `Rscript src/4.Datasets/data/simulated_annotations/generate_dataset.R`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'`
- `JWAS_SIMULATED_MT_SEEDS=11,22 JWAS_SIMULATED_MT_CHAIN_LENGTH=200 JWAS_SIMULATED_MT_BURNIN=50 JWAS_SIMULATED_MT_OUTPUT_FREQ=50 julia --project=. --startup-file=no benchmarks/simulated_annotations_multitrait_comparison.jl /tmp/jwas_mt_sim_smoke`
- `JWAS_SIMULATED_MT_SEEDS=101,202 JWAS_SIMULATED_MT_CHAIN_LENGTH=2000 JWAS_SIMULATED_MT_BURNIN=500 JWAS_SIMULATED_MT_OUTPUT_FREQ=20 julia --project=. --startup-file=no benchmarks/simulated_annotations_multitrait_comparison.jl /tmp/jwas_mt_sim_benchmark_20260407`
- `julia --project=docs --startup-file=no docs/make.jl`

## Result

The repository now includes a reproducible 2-trait annotation benchmark fixture
and a production JWAS benchmark script that can be rerun to compare the new
dense annotated 2-trait BayesC method against the requested baselines.
