# BayesR Fast Blocks Long-Chain Schedule Comparison Implementation

## What Changed

Added a new benchmark mode to `benchmarks/bayesr_fast_blocks_parity.jl`:

- `JWAS_BAYESR_BLOCK_MODE=long_chain_schedule_comparison`

The mode runs three BayesR schedules over multiple seeds:

1. production dense BayesR
2. benchmark-local default-blocks single-rep BayesR
3. production burnin-gated BayesR fast blocks

It writes:

- `schedule_runs.csv`
- `schedule_pairwise_summary.csv`

Also added a subprocess regression test to `test/unit/test_bayesr_parity.jl` that exercises the new mode with a small multiseed configuration.

## Files Changed

- `benchmarks/bayesr_fast_blocks_parity.jl`
- `test/unit/test_bayesr_parity.jl`
- `docs/plans/2026-03-20-bayesr-fast-blocks-long-chain-schedule-comparison-design.md`
- `docs/plans/2026-03-20-bayesr-fast-blocks-long-chain-schedule-comparison-plan.md`

## Verification

Ran:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'`
- `JWAS_BAYESR_BLOCK_MODE=long_chain_schedule_comparison JWAS_BAYESR_BLOCK_CHAIN_LENGTH=10000 JWAS_BAYESR_BLOCK_BURNIN=2000 JWAS_BAYESR_BLOCK_SEEDS=2026,2027,2028,2029,2030 julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl /tmp/bayesr_fast_blocks_long_chain_schedule_comparison`
- `julia --project=. --startup-file=no test/runtests.jl`

## Results

The long-chain benchmark gave two different answers:

### 1. Production burnin-gated fast blocks versus dense

This comparison was reasonably close over `10000/2000` and five seeds:

- mean `sigmaSq` abs diff: `0.0435`
- max `sigmaSq` abs diff: `0.0914`
- mean residual-variance abs diff: `0.00431`
- mean nonzero-frequency abs diff: `0.0362`
- mean `pi_class1` abs diff: `0.0334`

So the long chain does reduce the dense-vs-fast-block gap to a fairly small level.

### 2. Benchmark-local single-rep default blocks versus production dense

This comparison did **not** stay close over long chains:

- mean `sigmaSq` abs diff: `1.02`
- mean nonzero-frequency abs diff: `0.106`
- large `pi` shifts

This is a benchmark caveat, not evidence that `nreps = 1` is bad. The benchmark-local single-rep Gibbs loop is not a production-equivalent oracle over long chains. The earlier exact single-rep match came from comparing two local loops, not local single-rep against production dense.

## Interpretation

The main actionable result is:

- production burnin-gated BayesR fast blocks are much closer to dense BayesR on long chains than they looked on shorter diagnostics

The remaining open question is not the long-chain burnin-gated schedule itself. It is whether we want a production-quality way to test `default blocks + nreps = 1` against production dense without relying on the current benchmark-local Gibbs driver.
