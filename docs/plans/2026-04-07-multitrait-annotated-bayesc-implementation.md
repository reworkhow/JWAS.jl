# Multi-Trait Annotated BayesC Implementation

## Summary

Implemented dense 2-trait annotated BayesC on the production JWAS path using a
3-step tree annotation prior:

- step 1: `00` versus active
- step 2: `11` versus singleton `{10, 01}`
- step 3: `10` versus `01`

The implementation keeps the existing dense multi-trait BayesC likelihood and
marker-effect updates, adds marker-specific 4-state annotation priors for the
2-trait case, and preserves current single-trait annotated BayesC behavior.

## Code Changes

### Annotation setup and validation

Updated [readgenotypes.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/1.JWAS/src/markers/readgenotypes.jl),
[build_MME.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/1.JWAS/src/build_MME.jl),
and [input_data_validation.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/1.JWAS/src/input_data_validation.jl):

- delayed annotated BayesC `Pi` specialization so multi-trait startup can retain
  a joint prior dictionary
- added 2-trait annotated BayesC annotation-state initialization with
  `nsteps = 3`, `nclasses = 4`, zero coefficients, and repeated startup
  `snp_pi[m, 4]`
- added explicit v1 guardrails for unsupported combinations:
  `ntraits != 2`, `storage != :dense`, `fast_blocks = true`,
  `constraint = true`, and `RRM = true`

### Annotation Gibbs updates

Updated [annotation_updates.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/1.JWAS/src/MCMC/annotation_updates.jl):

- added 2-trait BayesC tree-step indicator helpers for:
  - `delta != 00`
  - `delta == 11` within active markers
  - `delta == 10` within singleton markers
- reused the existing binary probit annotation Gibbs machinery step-by-step
- rebuilt the per-marker 4-state prior matrix after each annotation update
- refreshed the user-facing `Mi.π` summary as the average marker-level joint
  prior over `00`, `10`, `01`, and `11`

### Multi-trait dense sampler integration

Updated [MCMC_BayesianAlphabet.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl)
and [MTBayesABC.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl):

- removed the old single-trait-only annotation stop for the supported 2-trait
  dense BayesC path
- routed annotated 2-trait BayesC through marker-specific prior lookups from
  `annotations.snp_pi`
- added annotated variants of sampler I and sampler II that replace shared
  `BigPi` lookups with per-marker 4-state priors while leaving the Gaussian
  likelihood calculations intact
- updated the MCMC order so multi-trait annotated BayesC refreshes annotation
  priors after the marker sweep instead of sampling a shared `Pi`

### Startup variance helpers and output

Updated [tools4genotypes.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/1.JWAS/src/markers/tools4genotypes.jl)
and [output.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/src/1.JWAS/src/output.jl):

- added `genetic2marker` support for annotation-backed 2-trait 4-state startup
  priors
- labeled the 2-trait annotation coefficient output by tree step:
  `zero_vs_active`, `11_vs_singleton`, and `10_vs_01`

## Tests

Updated [test_annotated_bayesc.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/test/unit/test_annotated_bayesc.jl):

- verified 2-trait annotated BayesC builds a 3-step annotation state
- verified the startup `snp_pi` rows match the supplied joint prior
- verified unsupported multi-trait paths fail with explicit errors

Updated [test_multitrait_mcmc.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/test/unit/test_multitrait_mcmc.jl):

- added an end-to-end dense 2-trait annotated BayesC run
- verified annotation coefficient and `pi_<geno>` outputs are produced

Updated the manuals:

- [annotated_bayesc.md](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/docs/src/manual/annotated_bayesc.md)
- [bayesc_bayesr_comparison.md](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/docs/src/manual/bayesc_bayesr_comparison.md)
- [workflow.md](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/docs/src/manual/workflow.md)

## Verification

Run:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_advanced_coverage.jl")'`
- `julia --project=. --startup-file=no test/runtests.jl`
- `julia --project=docs --startup-file=no docs/make.jl`

## Result

JWAS now supports dense 2-trait annotated BayesC on the production path with a
marker-specific 4-state annotation prior and explicit unsupported-path checks.
Single-trait annotated BayesC remains unchanged.
