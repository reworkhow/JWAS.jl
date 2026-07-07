# BayesR Class Frequency Output Implementation

## Summary

This change keeps the existing BayesR `Model_Frequency` behavior unchanged:
it is still the posterior nonzero marker frequency, `Pr(delta > 1)`.

BayesR marker-effect output now also includes exact class-specific posterior
frequencies:

- `Class1_Frequency = Pr(delta == 1)`
- `Class2_Frequency = Pr(delta == 2)`
- `Class3_Frequency = Pr(delta == 3)`
- `Class4_Frequency = Pr(delta == 4)`
- `MediumLarge_Frequency = Pr(delta >= 3)`

This lets users distinguish high overall BayesR PIP caused by class-2
small-effect assignments from PIP driven by medium or large classes.

## MCMC Sample Output

`runMCMC` now accepts `output_marker_effect_samples::Bool=true`.

When `output_marker_effect_samples=false`, JWAS skips the large
`MCMC_samples_marker_effects_*` files. It still computes final marker-effect
summaries, including the new BayesR class-frequency columns, and still writes
smaller sample files such as:

- `MCMC_samples_marker_effects_variances_*`
- `MCMC_samples_pi_*`
- residual variance samples
- EBV and heritability samples when those outputs are enabled

This option is rejected with `causal_structure`, because SEM post-processing
requires the direct marker-effect MCMC sample files to generate indirect and
overall marker-effect samples.

The implementation deliberately does not use `output_samples_frequency=0` as
the public off switch, because that setting is currently tied to posterior
summary accumulation cadence.

## Files Changed

- `src/1.JWAS/src/JWAS.jl`
- `src/1.JWAS/src/types.jl`
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- `src/1.JWAS/src/output.jl`
- `test/unit/test_bayesr.jl`
- `test/unit/test_sem_comprehensive.jl`
- `docs/src/manual/bayesc_bayesr_comparison.md`
- `docs/src/manual/workflow.md`

## Verification

Focused BayesR test:

```bash
/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Result: PASS.

The focused test covers:

- default BayesR marker-effect sample files remain written
- `output_marker_effect_samples=false` suppresses `MCMC_samples_marker_effects_*`
- marker-effect variance and pi sample files are still written
- BayesR marker output contains all class-frequency columns
- `Model_Frequency == Class2_Frequency + Class3_Frequency + Class4_Frequency`
- `MediumLarge_Frequency == Class3_Frequency + Class4_Frequency`

SEM validation test:

```bash
/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia --project=. --startup-file=no -e 'include("test/unit/test_sem_comprehensive.jl")'
```

Result: PASS after adding the explicit validation that
`output_marker_effect_samples=false` is not supported with `causal_structure`.

Annotated BayesR regression test:

```bash
/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Result: PASS.

Full test suite:

```bash
/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia --project=. --startup-file=no test/runtests.jl
```

Result: PASS, 988 tests in 2m40.1s.

Docs build:

```bash
/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia --project=docs --startup-file=no docs/make.jl
```

Result: blocked in this local environment. The first run failed because the
docs environment was not instantiated. After instantiating, the CI-style local
package setup with `Pkg.develop(PackageSpec(path=pwd()))` could not resolve
under Julia 1.9.3 because local JWAS dependency compatibility requires a newer
`Revise` range than Julia 1.9 can satisfy. No Julia 1.11 binary was available
locally.

## Notes

Raw per-iteration BayesR delta class samples are still not saved. The added
class-frequency columns are exact posterior summaries accumulated during MCMC,
which avoids adding another large marker-by-sample output file by default.
