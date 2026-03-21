# BayesR Fast Blocks Long-Chain Schedule Comparison Design

## Goal

Run a long-chain benchmark that compares BayesR fast-block schedules against dense BayesR on the same synthetic dataset across multiple seeds.

## Question

The open question is whether the remaining BayesR fast-block drift is mainly:

- a short-chain stabilization issue, or
- a persistent effect of repeated within-block sweeps

We already know:

- `block_size = 1` matches dense BayesR
- default block partition plus `nreps = 1` for all iterations also matches dense BayesR

So the relevant long-chain comparison is now the schedule itself.

## Scope

Benchmark-only change.

Compare these three methods:

1. dense BayesR
2. default block partition with `nreps = 1` for all iterations
3. production BayesR fast-block schedule
   - `nreps = 1` during burnin
   - `nreps = block_size` after burnin

Use:

- same synthetic dataset
- `chain_length = 10000`
- `burnin = 2000`
- seeds `2026,2027,2028,2029,2030`
- free `pi`
- free `sigmaSq`

## Recommendation

Use a new mode in the existing `benchmarks/bayesr_fast_blocks_parity.jl` harness rather than an ad hoc one-off script. That keeps the experiment reproducible and lets later schedule changes reuse the same comparison.

## Outputs

Write:

- `schedule_runs.csv`
- `schedule_pairwise_summary.csv`

The per-seed run table should include:

- method
- seed
- sigmaSq
- residual variance
- mean nonzero frequency
- pi class values
- runtime

The pairwise summary should compare:

- `single_rep` vs `dense`
- `burnin_gated` vs `dense`

for mean and max differences across seeds.

## Decision Rule

- If `burnin_gated` stays close to `single_rep` and `dense` over long chains, then the earlier drift was largely a short-chain stabilization issue.
- If `burnin_gated` still differs materially while `single_rep` remains close to dense, then the repeated within-block sweeps are the source of the remaining BayesR fast-block discrepancy.
