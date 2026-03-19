# BayesR Controlled Replay Report

## Setup

- Benchmark: one-iteration controlled replay
- Mode: fixed `pi`
- Initial state: exported shared state from `benchmarks/bayesr_parity_common.jl`
- Shared draws: `data/replay_draws_iteration1.csv`
- JWAS runner: `benchmarks/bayesr_parity_replay_jwas.jl`
- R runner: `benchmarks/bayesr_parity_replay_reference.R`
- Comparison: `benchmarks/bayesr_parity_replay_compare.jl`

## Result

The controlled replay does **not** show a BayesR logic mismatch at iteration 1.

With the same exported draws:

- `chosen_class` matched for every marker
- `new_alpha` matched for every marker
- `ssq` matched exactly
- `nnz` matched exactly
- `sigmaSq_new` matched to numerical precision (`2.27e-12` abs diff)

The remaining differences are small deterministic numeric gaps:

- first scalar mismatch: `mu_old`
- first marker mismatch: `rhs` at `m1`
- max scalar abs diff: `8.819169e-06` on `vare_new`
- max marker abs diff: `5.087024e-06`

## Interpretation

The replay result shifts the parity question:

- the BayesR update formulas are aligned across JWAS and the R reference for iteration 1 when both consume the same draws
- the previous raw-trace divergence is not caused by an immediate formula mismatch in class probabilities, class choice, marker-effect sampling, or `sigmaSq` sampling
- the visible replay differences are consistent with numeric precision differences in the JWAS benchmark-local path, which currently follows the genotype backend element type used by JWAS

So the current evidence says:

- there is no one-step BayesR implementation bug exposed by controlled replay
- the longer-run posterior `sigmaSq` discrepancy must come from later chain behavior, accumulated numerical precision effects, or another path outside this one-step replay

## Artifacts

- `/tmp/bayesr_controlled_replay_run/jwas_fixed_pi/replay_marker_iteration1.csv`
- `/tmp/bayesr_controlled_replay_run/jwas_fixed_pi/replay_scalar_iteration1.csv`
- `/tmp/bayesr_controlled_replay_run/ref_fixed_pi/replay_marker_iteration1.csv`
- `/tmp/bayesr_controlled_replay_run/ref_fixed_pi/replay_scalar_iteration1.csv`
- `/tmp/bayesr_controlled_replay_run/comparison_replay_marker_iteration1.csv`
- `/tmp/bayesr_controlled_replay_run/comparison_replay_scalar_iteration1.csv`
