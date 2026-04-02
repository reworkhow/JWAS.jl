# Annotated BayesR Jian-Style Startup Implementation

## Summary

Dense annotated BayesR startup now matches Jian's startup structure more closely:

- `snp_pi` starts as the supplied common class prior repeated across markers
- annotation coefficients start at zero
- no annotation pre-fit happens before the first phenotype-informed BayesR sweep

## Production Changes

### `src/1.JWAS/src/markers/readgenotypes.jl`

- kept the BayesR nondegenerate conditional-split validation
- removed the startup intercept initialization from the starting `Pi`
- annotated BayesR annotation objects now start with zero coefficients and the
  supplied constant `snp_pi`

### `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

- removed the pre-MCMC annotation initialization/update for annotated BayesR
- the first annotation update now happens only in the normal Gibbs loop after the
  first marker sweep

## Tests

Updated `test/unit/test_annotated_bayesr.jl` to check:

- annotated BayesR startup uses zero coefficients
- startup `mu` is zero
- startup `snp_pi` matches the supplied starting `Pi`
- startup is deterministic across RNG seeds

Regression checks also still pass for:

- annotated BayesC unit tests
- BayesR unit tests

## Result

The Jian-style startup materially improved dense annotated BayesR cross-seed
agreement on the large reproducer, but it did not fully solve the instability.
The remaining problem is therefore no longer just startup mismatch.
