# Simulated Annotations Cross-Validation Design

## Goal

Add an out-of-fold genomic prediction benchmark on the packaged
`simulated_annotations` fixture using the production JWAS path. The new
benchmark should compare:

- current multi-trait BayesC-family methods
- single-trait annotated BayesC
- single-trait annotated BayesR
- single-trait BayesC
- single-trait BayesR

The primary purpose is genomic prediction under cross-validation rather than
in-sample `cor(y, EBV)`.

## Recommended Design

Use individual-level K-fold cross-validation with the same fold assignment for
every method. In each held-out fold:

- multi-trait runs hide both `y1` and `y2` for held-out individuals
- single-trait runs hide the target trait for held-out individuals
- the full genotype matrix remains available for all individuals
- EBV output is requested for all individuals, but only the held-out fold is
  scored for out-of-fold prediction

This is the cleanest apples-to-apples comparison because the multi-trait models
do not get extra observed test-fold trait information that the single-trait
models cannot use.

## Approaches Considered

### 1. Standard K-fold CV with both traits hidden in the test fold

This is the recommended approach.

Pros:

- direct genomic prediction benchmark
- consistent across multi-trait and single-trait methods
- reuses the current production benchmark harness cleanly
- easy to summarize fold-wise and seed-wise

Cons:

- multi-trait models are not allowed to exploit cross-trait information in the
  test fold itself

### 2. Trait-masked CV for multi-trait methods

Example: keep `y1` observed in the test fold and predict `y2`, then reverse.

Pros:

- measures conditional prediction strength from one trait to another

Cons:

- not comparable to the requested single-trait baselines
- would require a second benchmark mode and separate reporting

### 3. Repeated random holdout instead of K-fold

Pros:

- simple to implement
- easy to parallelize later

Cons:

- more variable fold composition
- less standard than K-fold for this small packaged benchmark
- harder to compare fold-by-fold across methods

## Method Set

The CV benchmark should include:

- `MT_BayesC`
- `MT_BayesC_I`
- `MT_BayesC_II`
- `MT_Annotated_BayesC_I`
- `MT_Annotated_BayesC_II`
- `MT_EmptyAnnotated_BayesC_I`
- `MT_EmptyAnnotated_BayesC_II`
- `BayesC_single`
- `Annotated_BayesC_single`
- `BayesR_single`
- `Annotated_BayesR_single`

## Fold Structure

Recommended default:

- `5` folds
- deterministic fold assignment from individual IDs plus a benchmark seed
- same folds reused across all methods within a benchmark repetition

To preserve the existing seed-based benchmark pattern, each benchmark seed
should define one fold partition. Final summaries should average over:

- fold within seed
- seed across replicated partitions

## Data Flow

For each benchmark seed:

1. create fold assignments on individual IDs
2. for each fold:
   - write a masked phenotype table for the held-out fold
   - run the method through the production benchmark path
   - collect EBV output for all IDs
   - score only the held-out IDs
3. aggregate fold-level prediction metrics into seed-level summaries
4. aggregate seed-level summaries into method summaries

The current marker-level recovery and pleiotropy analyses should remain
available for the non-CV benchmark. The CV extension should focus on
out-of-fold prediction metrics.

## Prediction Metrics

Primary CV metrics:

- `cor(y1, EBV_y1)` on held-out individuals
- `cor(y2, EBV_y2)` on held-out individuals

Secondary CV metrics:

- RMSE on held-out individuals for each trait
- mean prediction correlation across traits

For single-trait methods, store one row per trait and also report a pooled
family-level mean across `y1` and `y2`.

## Runtime And Output

The benchmark should save:

- per-fold prediction tables joined by `ID`
- per-fold metrics
- per-seed method summaries
- final method summary

Suggested filenames:

- `cv_fold_assignments.csv`
- `cv_per_fold_summary.csv`
- `cv_method_summary.csv`

## Validation And Testing

Add coverage for:

- deterministic fold assignment
- no overlap between train and test IDs within a fold
- union of fold test IDs equals the full ID set
- EBV/phenotype joins preserve the held-out row count

Benchmark verification should include:

- targeted CV coverage test
- full `test/runtests.jl` after implementation

## Success Criteria

The change is complete when:

- the production benchmark can run K-fold CV on the packaged fixture
- all requested method families produce out-of-fold prediction summaries
- the report clearly separates CV prediction results from the earlier
  in-sample and marker-recovery benchmarks
