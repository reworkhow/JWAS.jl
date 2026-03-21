# BayesR Fast-Blocks Benchmark Cleanup Design

## Goal

Clean up the BayesR fast-blocks benchmark layer without changing the production
BayesR sampler.

## Design

### Production code

Do not change:

- `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- `src/1.JWAS/src/JWAS.jl`

The current production BayesR block behavior is accepted for now.

### Benchmark split

Keep the main benchmark script focused on production-facing comparisons:

- fixed-hyperparameter production comparison
- dense multiseed summary
- long-chain production schedule comparison

Move exploratory benchmark-local diagnostics out of the main script:

- local single-rep loops
- local trace loops
- exploratory debug-only modes

Those diagnostics will live in a separate debug script under `benchmarks/debug/`.

### Long-chain schedule mode

Simplify `long_chain_schedule_comparison` so it compares only:

- `dense`
- `burnin_gated`

The benchmark-local `single_rep` branch should not remain in the production-facing
summary mode because it is not a reliable production oracle over long chains.

### Tests

Keep subprocess coverage for the production-facing benchmark mode.

Remove or stop requiring subprocess coverage for debug-only modes from the main
benchmark test file.

### Reports

Keep the new summary reports in `benchmarks/reports/`.

These are durable investigation artifacts and should not be tied to whether the
debug script remains in the main benchmark path.
