# Annotated BayesR Benchmark Design

## Goal

Benchmark dense single-trait annotated BayesR on the production JWAS path and
compare it against ordinary dense BayesR and annotated BayesC on a synthetic
dataset with informative annotations.

## Scope

Included:

- production `runMCMC` path only
- single trait
- dense storage only
- methods:
  - `BayesR`
  - annotated `BayesR`
  - annotated `BayesC`
- long-chain multiseed comparison
- synthetic data where annotations are part of the data-generating process
- report with method-level and truth-aware summaries

Excluded:

- `fast_blocks`
- `storage=:stream`
- multi-trait
- real-data benchmark
- annotation-contribution-to-phenotype summaries

## Benchmark Question

The benchmark should answer:

- whether annotated BayesR improves SNP prioritization relative to ordinary
  BayesR when annotations truly enrich nonzero and larger-effect classes
- how annotated BayesR compares to annotated BayesC on posterior PIP and EBV
  quality under the same synthetic truth
- whether the fitted annotation coefficients recover the expected direction of
  enrichment

## Synthetic Data

Use one synthetic scenario with informative annotations.

Recommended defaults:

- `n_obs = 200`
- `n_markers = 1000`
- two annotation columns:
  - `annotation_1`: informative
  - `annotation_2`: nuisance or null
- moderate trait heritability around `0.45`

Data-generating idea:

- simulate genotype matrix with marker values in `0/1/2`
- assign each SNP a true latent mixture class
- let `annotation_1` increase the probability of:
  - being nonzero
  - being in the larger-effect BayesR classes
- let `annotation_2` have little or no effect
- sample marker effects from the class-specific BayesR variances
- generate phenotype by `y = X * beta + e`

This ensures that the annotation signal is real and benchmarkable.

## Methods

Run these production methods on the same dataset for each seed:

- dense `BayesR`
- dense annotated `BayesR`
- dense annotated `BayesC`

Only the presence of annotations and the method change. Genotypes, phenotype,
chain settings, and seeds stay fixed across methods.

## Chain Settings

Use long enough chains to compare posterior summaries rather than short-chain
noise.

Recommended defaults:

- `chain_length = 10000`
- `burnin = 2000`
- `output_samples_frequency = 10`
- seeds `2026:2030`

The smoke test can override these through environment variables with much
smaller settings.

## Metrics

Primary metrics:

- `cor(y, EBV)`
- posterior mean PIP
- posterior mean PIP among:
  - truly causal SNPs
  - null SNPs
  - annotation-enriched SNPs
  - non-enriched SNPs
- top-marker ranking against truth

Annotated BayesR-specific summaries:

- posterior mean class occupancy
- posterior mean `pi`
- annotation coefficient posterior means and SDs

Annotated BayesC-specific summaries:

- posterior mean inclusion probability
- annotation coefficient posterior means and SDs

## Outputs

The benchmark script should write:

- `comparison_runs.csv`
  - one row per method and seed
- `comparison_summary.csv`
  - pooled method summaries across seeds
- `pip_group_summary.csv`
  - grouped PIP summaries by truth label and annotation group
- `annotation_coefficients.csv`
  - posterior coefficient summaries for annotated methods
- `truth_metadata.csv`
  - per-marker truth, annotation labels, and true class

## Files

Create:

- `benchmarks/annotated_bayesr_comparison.jl`
- `benchmarks/reports/2026-03-28-annotated-bayesr-benchmark-report.md`
- `docs/plans/2026-03-28-annotated-bayesr-benchmark-implementation.md`

Modify:

- `test/unit/test_bayesr_parity.jl`

## Verification

Required verification:

- subprocess smoke test for the benchmark script
- targeted benchmark run on the full intended settings
- `julia --project=. --startup-file=no test/runtests.jl`

If the benchmark report is added to docs only as a file and not to the built
manual, `docs/make.jl` is not required for this task.
