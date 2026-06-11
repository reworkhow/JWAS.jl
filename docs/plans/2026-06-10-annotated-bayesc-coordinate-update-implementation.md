# Annotated BayesC Coordinate Update Implementation

## Code Change

Single-trait annotated BayesC now updates the binary annotation
coefficients with the existing coordinate Gibbs helper:

- intercept: flat prior
- annotation slopes: normal prior with variance `ann.variance`

The runtime path changed in
`src/1.JWAS/src/MCMC/annotation_updates.jl`:

```julia
latent_residual = ann.liability .- ann.mu
gibbs_update_binary_probit_annotation_coefficients!(
    ann.coefficients,
    ann.design_matrix,
    latent_residual,
    ann.variance,
)
ann.mu .= ann.design_matrix * ann.coefficients
```

No single-trait annotation-slope variance sampler was added. That remains a
separate model change.

## Tests

Added a unit test that verifies single-trait annotated BayesC uses the
coordinate probit update with a flat intercept and shrunken slopes. The test was
run before the production change and failed against the old joint Gibbs update:

- `annotation sampler uses standard probit latent variance`: coefficient and
  `Pi` mismatches
- `single-trait annotated BayesC uses coordinate probit update with shrunken
  slopes`: coefficient, `mu`, and `Pi` mismatches

After the production change:

```bash
julia --project=. --startup-file=no test/unit/test_annotated_bayesc.jl
```

passed:

- `Annotated BayesC API and validation`: 104/104
- `Standard BayesC preserves scalar pi output`: 4/4
- `Annotated BayesC dense run`: 6/6
- `Annotated BayesC fast_blocks run`: 4/4
- `Annotated BayesC independent fast_blocks run`: 6/6
- `Annotated BayesC streaming run`: 4/4

## Test Data Benchmark

Command:

```bash
julia --project=/tmp/jwas_hdf5_env --startup-file=no benchmarks/single_trait_test_data_coordinate_update_check.jl
```

Settings:

- data: `/Users/haocheng/Downloads/test_data`
- chain length: 2000
- burnin: 500
- output sample frequency: 20
- seed: 101
- starting `Pi`: 0.99
- SNPs: 49,224
- individuals: 500

Outputs:

`benchmarks/reports/single_trait_test_data_pip_comparison_after_coordinate_update/`

### PIP Summary

| Run | Mean PIP | Median PIP | Max PIP | Fraction PIP <= 0.01 | Fraction PIP > 0.1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| after init fix, full annotated | 0.006685 | 0.000000 | 0.160000 | 0.643832 | 0.000366 |
| after coordinate update, full annotated | 0.009420 | 0.013333 | 0.093333 | 0.497908 | 0.000000 |
| intercept-only annotated control | 0.010309 | 0.013333 | 0.080000 | 0.460568 | 0.000000 |

The intercept-only control is unchanged by the coordinate update, as expected,
because it has no slopes.

### Full Annotated Trace Summary

| Run | Samples | First | Middle | Last | Mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| after init fix, full annotated | 75 | 0.006314 | 0.006673 | 0.006878 | 0.006791 |
| after coordinate update, full annotated | 75 | 0.008189 | 0.010556 | 0.006761 | 0.009477 |
| intercept-only annotated control | 75 | 0.006922 | 0.012990 | 0.011806 | 0.010408 |

### Annotation Coefficients

After the initialization fix only, the short full-annotation run still had large
annotation slope excursions:

| Coefficient | Estimate | SD |
| --- | ---: | ---: |
| Intercept | -2.516981 | 0.073497 |
| Annotation_1 | -1.028352 | 0.632878 |
| Annotation_2 | 0.649881 | 0.394985 |
| Annotation_3 | -0.055161 | 0.420647 |
| Annotation_4 | -1.993577 | 0.845473 |

After the coordinate update:

| Coefficient | Estimate | SD |
| --- | ---: | ---: |
| Intercept | -2.331558 | 0.039720 |
| Annotation_1 | -1.036212 | 0.482253 |
| Annotation_2 | -0.057809 | 0.216092 |
| Annotation_3 | -0.181507 | 0.248847 |
| Annotation_4 | -0.012742 | 0.202317 |

The intercept is now close to the intercept-only control estimate
(`-2.318409`), while most annotation slopes are pulled much closer to zero.

## Interpretation

On this small test data, the coordinate update does not create PIP inflation.
It makes the full annotated short-chain result look more like the intercept-only
baseline overall, while reducing the large class-specific slope excursions seen
after the initialization-only fix.

This supports the diagnosis that the single-trait coefficient update was part
of the poor-mixing behavior for full annotated BayesC. The benchmark is still a
short-chain diagnostic, not a final convergence study.
