# Multi-Trait Annotated BayesR Implementation

## Summary

Implemented dense 2-trait annotated BayesR with a 16-state trait-by-magnitude prior, trait-specific annotation effects for BayesR magnitude classes, and a full global cross-trait marker-effect covariance matrix.

The implementation keeps the minimum-update design:

- reuse multi-trait annotated BayesC's 00/10/01/11 active-pattern tree
- add trait-specific nested BayesR magnitude steps
- keep latent `beta` sampled for every SNP and trait
- define realized marker effects as `alpha[t, j] = sqrt(gamma[class[t, j]]) * beta[t, j]`
- update `G` from all latent `beta` vectors across all SNPs

## Files Changed

- `src/1.JWAS/src/markers/annotation_setup.jl`
  - added 16-state BayesR helpers
  - added model-dependent 2-trait BayesR annotation initialization
- `src/1.JWAS/src/MCMC/annotation_updates.jl`
  - added 16-state prior rebuild
  - added seven nested probit response steps for multi-trait BayesR
- `src/1.JWAS/src/markers/tools4genotypes.jl`
  - initialized marker covariance from 16-state BayesR class priors
- `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesR.jl`
  - added dense 2-trait BayesR sampler with trait-wise class updates
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
  - dispatched dense multi-trait BayesR to `MTBayesR!`
- `src/1.JWAS/src/variance_components.jl`
  - sampled multi-trait BayesR `G` from all latent `beta` vectors
- `src/1.JWAS/src/output.jl`
  - added seven annotation step labels for multi-trait BayesR
- `test/unit/test_annotated_bayesr.jl`
  - added state, prior, response, covariance, and guardrail tests
- `test/unit/test_multitrait_mcmc.jl`
  - added conditional-helper and end-to-end dense 2-trait annotated BayesR tests
- `benchmarks/simulated_annotations_multitrait_comparison.jl`
  - added a focused `MT_Annotated_BayesR` simulated-annotations benchmark case
- `benchmarks/reports/2026-06-24-multitrait-annotated-bayesr-simulated-annotations-report.md`
  - recorded the focused simulated-annotations smoke run

## Verification

Focused tests run during implementation:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'`

Focused simulation smoke test:

- `JWAS_SIMULATED_MT_FOCUS_MODE=mt_annotated_bayesr JWAS_SIMULATED_MT_SEEDS=101 JWAS_SIMULATED_MT_CHAIN_LENGTH=300 JWAS_SIMULATED_MT_BURNIN=100 JWAS_SIMULATED_MT_OUTPUT_FREQ=20 JWAS_SIMULATED_MT_WARMUP=false julia --project=. --startup-file=no benchmarks/simulated_annotations_multitrait_comparison.jl benchmarks/out/mt_annotated_bayesr_smoke_20260624`

The implementation record is a plan document only. If manual docs are updated later, run:

- `julia --project=docs --startup-file=no docs/make.jl`

## Known Limitations

- dense 2-trait only
- no streaming support
- no random regression model support
- no constrained marker covariance support
- no BayesR fast-block sampler for multi-trait v1
- `G` is currently updated from all latent `beta` vectors, including SNPs in state 00
- future benchmark candidate: compare the current all-marker `G` update with an active-only update using SNPs where at least one trait has class `> 1`
