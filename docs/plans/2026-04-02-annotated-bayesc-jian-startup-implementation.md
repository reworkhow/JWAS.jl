# Annotated BayesC Jian-Style Startup Implementation

## Summary

Dense annotated BayesC startup now matches Jian's startup structure much more
closely:

- `Mi.π` starts as the supplied common prior expanded to marker level
- annotation coefficients start at zero
- no annotation pre-fit happens before the first phenotype-informed marker sweep

## Production Changes

### `src/1.JWAS/src/markers/readgenotypes.jl`

- added `annotation_starting_pi(...)` for annotated BayesC startup
- annotated BayesC now stores startup `Pi` as a marker-level exclusion vector
- removed the earlier intercept-based annotation initialization
- BayesC annotation objects now start with zero coefficients and zero `mu`

### `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

- pre-MCMC annotation initialization/update is now restricted to annotated BayesR
- annotated BayesC no longer fits annotation effects before the first marker sweep

### `src/1.JWAS/src/MCMC/annotation_updates.jl`

- BayesC annotation updates now overwrite the full marker-level `π_j` vector
- kept the deterministic BayesC indicator helper only as a fallback helper if it
  is invoked directly

### `src/1.JWAS/src/markers/tools4genotypes.jl`

- `genetic2marker` now supports BayesC marker-level `π_j` vectors

## Tests

Updated `test/unit/test_annotated_bayesc.jl` to check:

- annotated BayesC startup uses marker-level `π_j`
- startup coefficients and `mu` are zero
- startup is deterministic across RNG seeds
- standard-probit BayesC annotation update still works
- BayesC marker-level `π_j` is accepted by `genetic2marker`

## Result

This startup alignment materially improved dense annotated BayesC seed stability
on the large reproducer. The instability problem is no longer dominated by the
startup path.
