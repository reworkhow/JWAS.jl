# BayesR Parity Monte Carlo Follow-up

## Question

After the controlled replay showed no one-step BayesR formula mismatch, the remaining question was whether the earlier JWAS-vs-R `sigmaSq` gap reflected a real implementation issue or ordinary Monte Carlo noise from short chains.

## Findings

### 1. Production JWAS and benchmark-local Julia agree

On the same dataset, seed, fixed `pi`, `chain_length=1000`, and `burnin=200`:

- production JWAS `sigmaSq`: `6.329144`
- benchmark-local Julia `sigmaSq`: `6.327308`
- relative difference: about `0.03%`

So the remaining parity question is not caused by a gap between production JWAS BayesR and the benchmark-local Julia driver.

### 2. Float32 precision is not the explanation

With the same per-iteration draws reused inside Julia:

- `Float32` vs `Float64` posterior mean `sigmaSq` relative difference: `3.47e-8`
- `Float32` vs `Float64` posterior mean `vare` relative difference: `3.43e-7`

That is far too small to explain the earlier short-chain JWAS-vs-R gaps.

### 3. Short-chain single-run parity was noisy

For fixed `pi`, `chain_length=1000`, `burnin=200`, across seeds `2021:2030`:

- mean JWAS-vs-R `sigmaSq` relative difference: `9.93%`
- max JWAS-vs-R `sigmaSq` relative difference: `27.40%`
- mean `vare` relative difference: `1.24%`
- mean nonzero-frequency absolute difference: `0.00324`

But within-method spread was also large:

- JWAS `sigmaSq` mean: `7.8340`, SD: `2.0567`
- R `sigmaSq` mean: `8.7269`, SD: `2.8776`

So short single-chain parity was not a stable decision rule.

### 4. Longer chains shrink the gap materially

For fixed `pi`, `chain_length=10000`, `burnin=2000`:

- seed `2026`: `sigmaSq` relative difference `2.96%`
- seed `2030`: `sigmaSq` relative difference `0.62%`

For `estimate_pi`, `chain_length=10000`, `burnin=2000`:

- seed `2026`: `sigmaSq` relative difference `0.72%`
- seed `2030`: `sigmaSq` relative difference `1.86%`
- max `pi` absolute difference: about `0.014` to `0.016`

Residual variance remained close in these long runs:

- fixed `pi`: about `1.16%` to `1.20%`
- `estimate_pi`: about `0.85%` to `1.21%`

## Conclusion

The current evidence does not support a substantive BayesR implementation bug.

The earlier `sigmaSq` mismatch was largely a short-chain Monte Carlo effect. For parity work, single short chains such as `1000/200` are too noisy to be a reliable pass/fail benchmark.

## Recommended Benchmark Protocol

Use one of these:

- long single-chain parity:
  - `chain_length=10000`
  - `burnin=2000`
- multi-seed parity:
  - aggregate several seeds and compare summary distributions instead of one chain

The benchmark harness now includes [bayesr_parity_multiseed.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/bayesr-v1/benchmarks/bayesr_parity_multiseed.jl) to automate this workflow.
