# BayesR Parity Final Note

## Protocol

The final parity check uses the multiseed protocol from [bayesr_parity_multiseed.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/bayesr-v1/benchmarks/bayesr_parity_multiseed.jl):

- `chain_length=10000`
- `burnin=2000`
- seeds: `2026,2027,2028,2029,2030`
- modes:
  - `fixed_pi`
  - `estimate_pi`

Outputs were written under:

- `/tmp/bayesr_final_parity/fixed_pi`
- `/tmp/bayesr_final_parity/estimate_pi`

## Fixed `pi`

Multiseed summary:

- mean `sigmaSq` relative difference: `0.032145`
- max `sigmaSq` relative difference: `0.055008`
- mean residual-variance relative difference: `0.005274`
- max residual-variance relative difference: `0.012015`
- mean nonzero-frequency absolute difference: `0.001446`
- max nonzero-frequency absolute difference: `0.003063`

Sigma scale across seeds:

- JWAS mean `sigmaSq`: `8.518016`
- JWAS `sigmaSq` SD: `3.195692`
- R mean `sigmaSq`: `8.468974`
- R `sigmaSq` SD: `2.942238`

Interpretation:

- fixed-`pi` parity is good under the long-chain protocol
- the remaining cross-language gap is small relative to posterior scale and between-seed variation

## `estimate_pi`

Multiseed summary:

- mean `sigmaSq` relative difference: `0.015964`
- max `sigmaSq` relative difference: `0.043060`
- mean residual-variance relative difference: `0.006239`
- max residual-variance relative difference: `0.012098`
- mean nonzero-frequency absolute difference: `0.007581`
- max nonzero-frequency absolute difference: `0.013729`
- mean max-`pi` absolute difference: `0.011302`
- max max-`pi` absolute difference: `0.015945`

Sigma scale across seeds:

- JWAS mean `sigmaSq`: `4.603474`
- JWAS `sigmaSq` SD: `1.500619`
- R mean `sigmaSq`: `4.535853`
- R `sigmaSq` SD: `1.382623`

Interpretation:

- `estimate_pi` parity is also good under the long-chain protocol
- mixture-weight and sparsity summaries remain close across implementations

## Final Conclusion

For BayesR v1, the parity evidence supports this statement:

- JWAS BayesR is numerically consistent with the R reference at the summary level when parity is evaluated with a long-chain or multiseed protocol
- the earlier short-chain mismatch was largely Monte Carlo noise, not evidence of a substantive BayesR implementation bug

This is the note to keep with the BayesR v1 parity work.
