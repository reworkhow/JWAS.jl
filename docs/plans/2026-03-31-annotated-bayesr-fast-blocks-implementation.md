# Annotated BayesR Fast-Blocks Implementation

## Goal

Enable `fast_blocks` for dense annotated BayesR and verify the supported path
with tests, documentation, and a production benchmark.

Design reference:

- `docs/plans/2026-03-31-annotated-bayesr-fast-blocks-design.md`

## What Changed

### 1. Validation

Removed the explicit validation guard that rejected annotated BayesR with
`fast_blocks != false`.

File:

- `src/1.JWAS/src/input_data_validation.jl`

This was the only production guard blocking the feature. The BayesR block
kernel already accepted annotation-driven per-marker priors via
`genotypes.annotations.snp_pi`.

### 2. Unit coverage

Replaced the old failure expectation with a real annotated BayesR block run.

File:

- `test/unit/test_annotated_bayesr.jl`

The new coverage verifies that:

- annotated BayesR runs successfully with `fast_blocks=true`
- marker effects and annotation coefficients are produced
- `snp_pi` remains a valid four-class probability matrix after the run

### 3. Production benchmark extension

Extended the existing annotated BayesR production benchmark to include:

- `Annotated_BayesR_fast_blocks_default`
- `Annotated_BayesR_fast_blocks_1`

File:

- `benchmarks/annotated_bayesr_comparison.jl`

Benchmark harness changes:

- added optional benchmark cases driven by
  `JWAS_ANNOT_BENCH_INCLUDE_FAST_BLOCKS=true`
- preserved numeric `fast_blocks=1` instead of accidentally coercing it to
  `true`
- made the benchmark burnin block-aware so block runs still retain post-burnin
  samples after the outer-chain length is shortened by `fast_blocks`
- added summary metadata:
  - `fast_blocks_setting`
  - `block_size`
  - `outer_chain_length`

### 4. Benchmark smoke test

Added smoke coverage for the new block variants.

File:

- `test/unit/test_bayesr_parity.jl`

The smoke test now verifies:

- the benchmark script emits the two annotated fast-block variants
- `fast_blocks=1` remains numeric with `block_size = 1`
- default `fast_blocks=true` uses a larger block size

### 5. Manual

Updated the annotated BayesR manual to describe dense block support and to make
the approximation caveat explicit.

File:

- `docs/src/manual/annotated_bayesr.md`

## Benchmark Result

Production benchmark report:

- `benchmarks/reports/2026-03-31-annotated-bayesr-fast-blocks-report.md`

Main conclusion:

- annotated BayesR with `fast_blocks` is now supported
- it is not a dense-parity mode
- default `fast_blocks=true` is materially faster than dense annotated BayesR
  on the tested scenario, but it also shifts posterior summaries

## Verification

Targeted verification run during implementation:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'`

Smoke benchmark run:

- `JWAS_ANNOT_BENCH_CHAIN_LENGTH=60`
- `JWAS_ANNOT_BENCH_BURNIN=20`
- `JWAS_ANNOT_BENCH_N_OBS=30`
- `JWAS_ANNOT_BENCH_N_MARKERS=40`
- `JWAS_ANNOT_BENCH_SEEDS=2026,2027`
- `JWAS_ANNOT_BENCH_SCENARIO=stepwise_annotation_signal`
- `JWAS_ANNOT_BENCH_INCLUDE_FAST_BLOCKS=true`
- `julia --project=. --startup-file=no benchmarks/annotated_bayesr_comparison.jl /tmp/annotated_bayesr_fast_blocks_smoke`

Long-chain production benchmark run:

- `JWAS_ANNOT_BENCH_INCLUDE_FAST_BLOCKS=true`
- `JWAS_ANNOT_BENCH_SCENARIO=stepwise_annotation_signal`
- `julia --project=. --startup-file=no benchmarks/annotated_bayesr_comparison.jl /tmp/annotated_bayesr_fast_blocks_benchmark_20260331`

Further verification still needed before merge:

- docs build on the final tree
- full `test/runtests.jl` attempt on the final tree
