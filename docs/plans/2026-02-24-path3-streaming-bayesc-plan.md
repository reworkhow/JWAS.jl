# Path 3 MVP Implementation Plan: Streaming Genotypes for Original BayesC

**Goal:** Add an opt-in streaming genotype backend for large data while keeping dense loading as the default, primary, and unchanged path.

## Locked Requirements

1. Existing dense genotype loading remains default and unchanged.
2. Streaming backend is explicit opt-in only (`storage=:stream`).
3. No deprecation or behavior change for dense users.

## Architecture Decision

### Recommended: Clear Separation + Shared BayesC Math Kernel

- Keep dense loader/sampler path intact.
- Add a separate streaming backend module and codec.
- Reuse one BayesC marker-update kernel between dense non-block and streaming non-block paths.

### Pros

1. Lowest regression risk for dense path.
2. Clear review boundaries (backend/IO changes isolated).
3. Avoids math drift between dense and streaming implementations.
4. Easier long-term maintenance/debugging.

### Cons

1. Small abstraction overhead.
2. Slightly more upfront design work than inline branching.

### Alternatives (not selected)

1. Inline `if storage==...` throughout existing code:
   - Pro: fewer new files initially.
   - Con: quickly tangles code and hurts reviewability.
2. Duplicate BayesC math for dense and streaming:
   - Pro: local optimization freedom.
   - Con: high bugfix duplication and divergence risk.

## Public API Changes

1. Add `prepare_streaming_genotypes(...)`.
2. Extend `get_genotypes(...; storage::Symbol=:dense)`:
   - `:dense` default, unchanged behavior.
   - `:stream` packed backend path.
3. Add backend metadata in `Genotypes`:
   - `storage_mode::Symbol`
   - `stream_backend`

## MVP Scope

1. Original BayesC only (`fast_blocks=false`).
2. Single-trait only.
3. Unit residual weights only.
4. Float32 only (`double_precision=false`).
5. Exact genotype/phenotype ID match and order required.
6. `outputEBV` and `output_heritability` deferred/disabled for streaming MVP.

## File-Level Tasks

1. New backend/codec module:
   - `src/1.JWAS/src/markers/streaming_genotypes.jl`
2. Wire includes and export:
   - `src/1.JWAS/src/JWAS.jl`
3. Extend `Genotypes` metadata:
   - `src/1.JWAS/src/types.jl`
4. Add `storage` branch in genotype loading:
   - `src/1.JWAS/src/markers/readgenotypes.jl`
5. Add streaming BayesC wrapper with shared update kernel:
   - `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl`
   - `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
6. Add streaming validation gates:
   - `src/1.JWAS/src/input_data_validation.jl`
7. Extend memory guard to streaming working set:
   - `src/1.JWAS/src/markers/tools4genotypes.jl`
   - `src/1.JWAS/src/JWAS.jl`
8. Docs updates:
   - `docs/src/manual/large_genotype_data_streaming.md`
   - `docs/src/manual/workflow.md`

## Testing Plan

1. Dense non-regression tests remain mandatory.
2. Add streaming codec tests (encode/decode, missing+centering behavior).
3. Add streaming-vs-dense BayesC equivalence test (fixed seed, small dataset).
4. Add streaming constraint tests (fast blocks, multi-trait, non-unit weights, double precision).
5. Keep streaming tests in dedicated testsets for clear CI signals.

## Acceptance Criteria

1. Dense path behavior unchanged by default.
2. Streaming enabled only with explicit `storage=:stream`.
3. Streaming BayesC runs without materializing dense `N x P` matrix.
4. Memory guard reports streaming `O(N+P)` marker working-set estimate.
5. Docs clearly state dense is primary/default and streaming is additive experimental.

## Assumptions

1. Complete genomic data only in streaming MVP.
2. One genotype category in streaming MVP.
3. Fast-block support is post-MVP.
