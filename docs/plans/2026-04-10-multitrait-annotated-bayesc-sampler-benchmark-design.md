# Multi-Trait Annotated BayesC Sampler Benchmark Design

## Goal

Design a careful comparison on the packaged `simulated_annotations` genotype fixture to evaluate:

- multi-trait annotated BayesC with Gibbs sampler I
- multi-trait annotated BayesC with Gibbs sampler II
- conventional multi-trait BayesC
- relevant single-trait BayesC and BayesR baselines

The main question is whether the new multi-trait annotated BayesC implementation, and especially the sampler-I versus sampler-II choice, improves recovery of pleiotropic loci while maintaining competitive genomic prediction accuracy.

## Dataset

Use the packaged `simulated_annotations` fixture in `src/4.Datasets/data/simulated_annotations`.

Required files:

- `genotypes.csv`
- `phenotypes_mt.csv`
- `annotations_mt.csv`
- `truth_mt.csv`

This keeps the comparison on a reproducible, versioned dataset already bundled with JWAS.

## Method Matrix

Run the following methods:

### Multi-trait methods

- conventional multi-trait BayesC
- multi-trait annotated BayesC with `multi_trait_sampler=:I`
- multi-trait annotated BayesC with `multi_trait_sampler=:II`

### Single-trait methods

Trait `y1` and trait `y2` should each be run separately for:

- BayesC
- annotated BayesC
- BayesR
- annotated BayesR

## Primary Metric

Primary headline metric: pleiotropic recovery.

### Multi-trait methods

Use the true joint posterior shared-state probability:

- `P11_j = Pr(delta_j = (1,1) | data)`

For each run:

- rank markers by `P11`
- declare the top `k_shared_true` markers as shared, where `k_shared_true` is the number of true shared loci in `truth_mt.csv`
- report shared precision, recall, F1, and counts

This keeps the shared comparison threshold-free apart from a truth-sized top-k rule.

### Single-trait methods

Do not use `min(PIP_y1, PIP_y2)` as the main shared metric.

Instead:

- detect top `k_y1_true` markers from the trait-1 run using trait-1 PIP
- detect top `k_y2_true` markers from the trait-2 run using trait-2 PIP
- define the declared shared set as the overlap of the two detected sets
- compare that overlap against the true shared `11` loci

This gives an intuitive single-trait pleiotropy proxy based on overlapping detected causal variants.

## Secondary Metrics

Report the following as secondary metrics:

- genomic prediction accuracy for each trait
- any-active recovery
- trait-specific top-k recall
- marker-effect correlation where available
- runtime

Prediction accuracy should be reported as the correlation between phenotype and EBV for each trait.

## Comparison Rules

Follow the repository comparison rules strictly:

- join marker outputs to truth by explicit `marker_id`
- normalize marker IDs before joining
- assert row counts, key-set equality, and join counts before summarizing
- save joined tables to disk for inspection

Do not compare rows by position.

## Benchmark Execution

Use multiple seeds and nontrivial chain lengths.

Recommended workflow:

1. smoke benchmark with short chains to validate the matrix
2. production benchmark with longer chains and multiple seeds
3. write a report under `benchmarks/reports/`

## Outputs

The updated benchmark should write:

- per-run summary
- method summary
- single-trait family any-active summary
- multi-trait shared posterior summaries
- single-trait shared-overlap summaries
- combined pleiotropy summaries

The report should highlight:

- sampler I versus sampler II for annotated multi-trait BayesC
- multi-trait annotated BayesC versus conventional multi-trait BayesC
- multi-trait methods versus single-trait baselines on shared recovery
