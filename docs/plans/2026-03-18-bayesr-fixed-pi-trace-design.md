# BayesR Fixed-Pi Trace Design

## Goal

Build a benchmark-only diagnostic that compares JWAS BayesR and the R reference iteration by iteration under fixed `pi`, using the same explicit initial state, so we can locate where the shared variance mismatch enters.

## Recommended Approach

1. `Benchmark-only trace extension`
   Extend the current parity harness so both JWAS and R write per-iteration traces for `sigmaSq`, `ssq`, `nnz`, and `vare` under fixed `pi`.
2. `Core sampler instrumentation`
   Add debug logging to production MCMC code and mirror it in R.
3. `One-step replay harness`
   Compare one Gibbs iteration at a time from saved states.

Recommendation: use the benchmark-only trace extension first. It keeps the investigation reproducible and avoids contaminating the production sampler with debug-only logic.

## Scope

- Fixed `pi` only
- Benchmark-only first
- Short chain: `100` iterations
- `burnin = 0`
- Same explicit initial state in both runners
- No production BayesR behavior changes in this step

## Trace Harness

Both runners should consume one shared exported dataset and write aligned trace files.

The fixed-`pi` trace file schema should include:

- `iter`
- `sigmaSq`
- `ssq`
- `nnz`
- `vare`

The comparison workflow should align the two trace files by iteration and report:

- per-iteration differences
- first iteration where a material gap appears
- summary of mean/median/max absolute differences

## Explicit Initial State

Do not rely on each runner's default initialization.

Write the initial state into the benchmark dataset/config and make both runners use it explicitly:

- `beta0 = zeros(m)`
- `delta0 = ones(Int, m)` for the zero-effect class
- `mu0 = mean(y)`
- `ycorr0 = y - mu0`
- `sigmaSq0 = start_sigma_sq`
- `vare0 = start_vare`
- fixed `pi0 = start_pi`

This explicit initialization is only for the parity trace harness. It is not a change to normal JWAS BayesR initialization.

## Diagnostic Goals

With fixed `pi` and identical initial state, the trace comparison should tell us which quantity diverges first:

1. `ssq` diverges first
   Then the mismatch is upstream in marker updates, class assignments, or `yCorr` bookkeeping.
2. `nnz` diverges first while `ssq/nnz` stays similar
   Then the main difference is class occupancy.
3. `ssq` and `nnz` stay close but `sigmaSq` diverges
   Then the issue is concentrated in the variance draw or hyperparameter interpretation.

Because Julia and R use different RNG streams, iteration-by-iteration identity is not the target. The target is to locate the drift pattern and identify which sufficient statistic moves first.

## Files

Expected benchmark-layer changes:

- Extend `benchmarks/bayesr_parity_common.jl`
- Extend `benchmarks/bayesr_parity_jwas.jl`
- Extend `benchmarks/bayesr_parity_reference.R`
- Extend or supplement `benchmarks/bayesr_parity_compare.jl`

Expected outputs:

- `data/initial_state.csv`
- `data/initial_scalars.csv` or equivalent config entries
- `jwas_fixed_pi/trace_fixed_pi.csv`
- `ref_fixed_pi/trace_fixed_pi.csv`
- `comparison_trace_fixed_pi.csv`

## Verification

- Add unit tests for the trace-file schema and initial-state export
- Run the fixed-`pi` trace benchmark end to end
- Keep the existing parity summary benchmark working
- Rerun the full Julia test suite before claiming completion

## Next Step

Implement the fixed-`pi` trace benchmark as a narrow diagnostic layer before touching the BayesR sampler again.
