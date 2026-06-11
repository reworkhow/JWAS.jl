# Single-Trait Annotated BayesC Annotation Variance Implementation

## Code Change

Single-trait annotated BayesC now samples scalar `ann.variance` after the
annotation probit coefficients are sampled, using the same scaled
inverse-chi-square update as Jian's `sbayesrc.R` and the existing JWAS nested
annotation sampler:

```text
sigma_alpha^2 = (sum(alpha_slopes^2) + 2) / chisq(number_of_slopes + 2)
```

The intercept is excluded from the sum of squares. Intercept-only annotation
models leave `ann.variance` unchanged.

The latent probit residual variance remains fixed at `1`; this change only
samples the prior variance for annotation slopes.

## Implementation

Updated `src/1.JWAS/src/MCMC/annotation_updates.jl`:

- factored the shared variance formula into
  `sample_annotation_effect_variance(coeffs)`
- kept the existing vector-per-step updater for nested annotation models
- added `sample_bayesc_binary_annotation_effect_variance!(ann)` for the
  scalar single-trait BayesC case
- called the scalar variance update after
  `gibbs_update_bayesc_binary_annotation_coefficients!(ann)` and before rebuilding
  marker-specific `Pi`

## Tests

Added unit coverage in `test/unit/test_annotated_bayesc.jl`:

- the single-trait annotated BayesC coordinate-update test now also checks the
  scalar variance draw from the post-update slopes
- a separate intercept-only test confirms scalar `ann.variance` is unchanged
  when there are no annotation slopes

Red test before implementation:

```text
actual_variance == expected_variance
0.25 == 0.3862911863063604
```

Green test after implementation:

```bash
julia --project=. --startup-file=no test/unit/test_annotated_bayesc.jl
```

passed:

- `Annotated BayesC API and validation`: 106/106
- `Standard BayesC preserves scalar pi output`: 4/4
- `Annotated BayesC dense run`: 6/6
- `Annotated BayesC fast_blocks run`: 4/4
- `Annotated BayesC independent fast_blocks run`: 6/6
- `Annotated BayesC streaming run`: 4/4
