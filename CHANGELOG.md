# Changelog

All notable changes to this project are documented in this file.

## 2026-02-26

### Added

- Low-memory, out-of-core streaming conversion for genotype text files in `prepare_streaming_genotypes`, including staged spool + transpose flow and disk guard checks.
- `conversion_mode` support for `:lowmem`, `:dense`, and `:auto` conversion paths.
- `auto_dense_max_bytes` threshold for automatic conversion-path selection.
- Unit coverage for low-memory conversion behavior, including dense/lowmem auto selection parity.

### Documentation

- Added usage guidance and examples for `conversion_mode` and `auto_dense_max_bytes`.
- Updated streaming backend and large-data walkthrough docs to reflect the low-memory conversion flow.
