# BayesR Fast Blocks Default-Blocks Single-Rep Implementation

## What Changed

Implemented a benchmark-only diagnostic mode in `benchmarks/bayesr_fast_blocks_parity.jl`:

- `JWAS_BAYESR_BLOCK_MODE=default_blocks_single_rep`

The mode:

- uses the default BayesR block partition
- keeps free `pi`
- keeps free `sigmaSq`
- keeps the same outer `chain_length` and `burnin` as dense BayesR
- forces one within-block sweep per outer iteration by calling `BayesR_block!` with a very large effective burnin

Also added a subprocess regression test in `test/unit/test_bayesr_parity.jl` that runs the benchmark mode and checks that the expected comparison CSVs are produced.

## Files Changed

- `benchmarks/bayesr_fast_blocks_parity.jl`
- `test/unit/test_bayesr_parity.jl`
- `docs/plans/2026-03-19-bayesr-fast-blocks-default-blocks-single-rep-design.md`
- `docs/plans/2026-03-19-bayesr-fast-blocks-default-blocks-single-rep-plan.md`

## Verification

Ran:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'`
- `JWAS_BAYESR_BLOCK_MODE=default_blocks_single_rep julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl /tmp/bayesr_fast_blocks_default_single_rep_run`

## Diagnostic Result

For the default benchmark dataset and seed, default block partition plus one within-block sweep matches dense BayesR essentially exactly:

- `sigmaSq` abs diff: `0.0`
- residual variance abs diff: `8.42e-7`
- mean nonzero frequency abs diff: `0.0`
- `pi` abs diff by class: `0.0`
- top marker-effect abs diffs: around `1e-7`

Interpretation:

- the default block partition itself is not the source of the remaining BayesR fast-block drift
- the remaining mismatch enters when repeated within-block sweeps are used
