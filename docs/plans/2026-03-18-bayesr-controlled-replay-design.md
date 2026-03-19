# BayesR Controlled Replay Design

## Goal

Build a benchmark-only one-iteration replay harness that forces JWAS BayesR and the R reference to consume the same random draws, so we can localize the first mismatch in the BayesR update path.

## Recommended Approach

1. `Benchmark-only controlled replay`
   Export one explicit draw list for iteration 1, then have Julia and R replay the same iteration from the same state and compare substeps.
2. `Debug helper in JWAS`
   Add production debug helpers that accept supplied draws.
3. `Hybrid benchmark replay`
   Keep the replay benchmark-only, but duplicate the BayesR update math there so every substep is visible.

Recommendation: use the hybrid benchmark replay first. It keeps production code untouched while making every deterministic and stochastic substep inspectable.

## Scope

- iteration 1 only
- fixed `pi` only
- benchmark-only
- same exported dataset
- same exported initial state
- same exported draw list
- no production BayesR behavior changes in this step

## Controlled Draws

Generate one shared draw file in Julia and export it to CSV.

The draw list should include:

- one normal draw for the `mu` update
- one uniform draw per marker for class assignment
- one normal draw per marker for the marker-effect draw
- one chi-square draw for `sigmaSq`
- one chi-square draw for `vare`

For class assignment, both runners should:

- compute the posterior class probabilities
- map the same `u_class ~ U(0,1)` through the cumulative probability vector

This avoids differences between Julia and R categorical samplers.

## Replay Outputs

Write two outputs per runner.

### 1. Marker Replay Table

One row per marker:

- `marker_id`
- `rhs`
- `old_alpha`
- `p_class1`
- `p_class2`
- `p_class3`
- `p_class4`
- `u_class`
- `chosen_class`
- `beta_hat_chosen`
- `inv_lhs_chosen`
- `z_beta`
- `new_alpha`
- `ycorr_norm_after`

### 2. Scalar Replay Table

One short scalar table:

- `mu_old`
- `mu_hat`
- `z_mu`
- `mu_new`
- `sigmaSq_old`
- `ssq`
- `nnz`
- `chisq_sigma`
- `sigmaSq_new`
- `vare_old`
- `chisq_vare`
- `vare_new`

### 3. Comparison Output

Align JWAS and R by marker ID and field name, then compute absolute differences.

## Comparison Order

The replay should be interpreted in this order:

1. `mu` update
2. per-marker deterministic quantities:
   - `rhs`
   - `inv_lhs`
   - `beta_hat`
   - class probabilities
3. per-marker stochastic application:
   - same `u_class` should give the same class
   - same `z_beta` should give the same `new_alpha`
4. end-of-sweep variance inputs:
   - `ssq`
   - `nnz`
5. final variance draws:
   - same `chisq_sigma` should give the same `sigmaSq_new`
   - same `chisq_vare` should give the same `vare_new`

If the first mismatch appears before any supplied draw is applied, the issue is deterministic math. If everything matches under shared draws for iteration 1, then the earlier divergence was driven by unsynchronized randomness rather than a one-step formula difference.

## Files

Expected benchmark-layer changes:

- extend `benchmarks/bayesr_parity_common.jl`
- add a JWAS replay runner under `benchmarks/`
- add an R replay runner under `benchmarks/`
- add a replay comparison script under `benchmarks/`
- extend parity unit tests for replay schema and comparison

Expected artifacts:

- `data/replay_draws_iteration1.csv`
- `jwas_fixed_pi/replay_marker_iteration1.csv`
- `jwas_fixed_pi/replay_scalar_iteration1.csv`
- `ref_fixed_pi/replay_marker_iteration1.csv`
- `ref_fixed_pi/replay_scalar_iteration1.csv`
- `comparison_replay_iteration1.csv`

## Verification

- add unit tests for draw export and replay comparison schema
- run the one-iteration replay end to end
- inspect the first mismatch location from the comparison output
- rerun targeted parity tests
- rerun the full Julia suite before claiming completion

## Expected Outcome

Either:

- we find a deterministic BayesR mismatch in a specific substep, or
- we show the two implementations agree exactly under shared draws for iteration 1

In either case, this turns the current broad parity question into a much narrower debugging target.
