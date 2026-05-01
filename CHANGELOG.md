# Changelog

All notable changes to this project are documented in this file.

## 2026-05-01

### Added

- `independent_blocks=true` support for fast-block marker samplers, with the exact sequential fast-block sweep remaining the default (`independent_blocks=false`).
- Explicit marker block starts for `fast_blocks`, including full-sweep chain semantics for user-supplied block partitions.
- Independent-block sampler coverage across single-trait BayesA/B/C, BayesR, annotated BayesC/BayesR, plain multi-trait BayesC, and dense 2-trait annotated BayesC paths.
- Packaged `simulated_annotations` dataset files for annotation-aware examples, including generated genotypes, phenotypes, annotations, truth data, raw genotype input, and the regeneration script.

### Documentation

- Added a dedicated dense 2-trait annotated BayesC manual page covering the 4-state prior, 3-step annotation parameterization, `Pi` validation, sampler selection, output interpretation, and current support restrictions.
- Refreshed annotated BayesC/BayesR and streaming genotype guidance to distinguish dense, fast-block, and streaming workflows.
- Documented the support boundary between exact sequential `fast_blocks` sweeps and approximate `independent_blocks=true` block-level parallelism.

## 2026-02-26

### Added

- Low-memory, out-of-core streaming conversion for genotype text files in `prepare_streaming_genotypes`, including staged spool + transpose flow and disk guard checks.
- `conversion_mode` support for `:lowmem`, `:dense`, and `:auto` conversion paths.
- `auto_dense_max_bytes` threshold for automatic conversion-path selection.
- Unit coverage for low-memory conversion behavior, including dense/lowmem auto selection parity.

### Documentation

- Added usage guidance and examples for `conversion_mode` and `auto_dense_max_bytes`.
- Updated streaming backend and large-data walkthrough docs to reflect the low-memory conversion flow.
