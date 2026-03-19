# BayesR Fixed-Pi Trace Report

## Setup

Commands used:

```bash
OUT=/tmp/bayesr_fixed_pi_trace_run
rm -rf "$OUT"
JWAS_PARITY_OUTDIR="$OUT" JWAS_PARITY_MODE=fixed_pi JWAS_PARITY_TRACE=1 JWAS_PARITY_CHAIN_LENGTH=100 JWAS_PARITY_BURNIN=0 julia --project=. --startup-file=no benchmarks/bayesr_parity_jwas.jl
Rscript benchmarks/bayesr_parity_reference.R "$OUT" fixed_pi trace
julia --project=. --startup-file=no benchmarks/bayesr_parity_compare.jl "$OUT" fixed_pi trace
```

Shared conditions:

- fixed `pi = [0.95, 0.03, 0.015, 0.005]`
- same exported dataset
- same explicit initial state
  - `beta0 = 0`
  - `delta0 = 1`
  - `mu0 = mean(y)`
  - `sigmaSq0 = start_sigma_sq`
  - `vare0 = start_vare`
- `chain_length = 100`
- `burnin = 0`

Compared trace columns:

- `sigmaSq`
- `ssq`
- `nnz`
- `vare`

## Main Findings

The raw iteration-by-iteration traces diverge immediately:

- first `sigmaSq` gap: iteration `1`
- first `vare` gap: iteration `1`
- first `ssq` gap: iteration `2`
- first `nnz` gap: iteration `2`

Aggregate differences over 100 iterations:

- mean `sigmaSq` abs diff: `10.3082`
- median `sigmaSq` abs diff: `4.8606`
- max `sigmaSq` abs diff: `118.3409`
- mean `ssq` abs diff: `14.3273`
- mean `nnz` abs diff: `1.06`

First rows of the comparison trace:

```text
iter 1: ssq diff = 0, nnz diff = 0, sigmaSq diff = 9.8571, vare diff = 0.1109
iter 2: ssq diff = 6.2650, nnz diff = 1, sigmaSq diff = 2.2456, vare diff = 0.3823
iter 3: ssq diff = 6.6005, nnz diff = 1, sigmaSq diff = 2.5359, vare diff = 0.1420
```

## Interpretation

This trace is informative, but not in the way originally hoped.

What it shows:

- With the same explicit initial state, both chains still separate immediately.
- At iteration `1`, `ssq` and `nnz` are still identical.
- The first mismatch is already in `sigmaSq` and `vare`.

What that means:

- The immediate iteration-1 gap is not evidence of a logic mismatch by itself.
- It is expected because Julia and R use different RNG streams for:
  - the mean update
  - the `sigmaSq` draw
  - the residual variance draw
- Once `sigmaSq` and `vare` differ at iteration `1`, the downstream marker updates see different conditional states, so `ssq` and `nnz` start diverging at iteration `2`.

So the fixed-`pi` raw trace comparison does **not** isolate the root cause of the posterior mean `sigmaSq` mismatch. It mostly confirms that unsynchronized random draws cause the two chains to decorrelate almost immediately.

## What This Rules Out

It does **not** support the claim that JWAS and R have a deterministic state-update mismatch at iteration `1`, because:

- `ssq` matches at iteration `1`
- `nnz` matches at iteration `1`
- both chains started from the same exported state

The first-step difference is already stochastic, not clearly structural.

## Recommended Next Step

Use a conditional or semi-deterministic diagnostic instead of raw trace equality.

Best next options:

1. Compare deterministic quantities given the same state:
   - class posterior probabilities
   - posterior means `betaHat`
   - `lhs` / `rhs`
2. Force selected random draws to be shared externally for one-step replay:
   - one `mu` draw
   - one `sigmaSq` chi-square draw
   - one `vare` chi-square draw
3. Compare trace distributions rather than paired iterations:
   - medians
   - quantiles
   - conditional summaries of `sigmaSq | nnz`

The most direct next debugging step is option `1`: compare deterministic BayesR update quantities for one iteration from the same saved state.
