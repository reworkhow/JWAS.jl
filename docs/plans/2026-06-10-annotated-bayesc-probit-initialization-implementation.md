# Annotated BayesC Probit Initialization Implementation

## Changes

- Added `initialize_bayesc_single_trait_annotation_coefficients!` in `src/1.JWAS/src/markers/annotation_setup.jl`.
- Single-trait annotated BayesC now initializes the annotation probit intercept from the starting marker inclusion probability `mean(1 - Pi)`.
- Annotation slopes still start at zero.
- Marker-level `genotypei.π` still stores the requested starting exclusion probabilities.
- Updated the annotated BayesC startup unit test to exercise the same finalization path used by `build_model`.

## Verification

Focused unit test:

```bash
julia --project=. --startup-file=no test/unit/test_annotated_bayesc.jl
```

Result: `100` tests passed in the annotated BayesC API and validation test set, plus the existing annotated BayesC smoke tests passed.

## Test-Data Impact Check

Data: `/Users/haocheng/Downloads/test_data`

Settings:

- chain length: `2000`
- burnin: `500`
- output samples frequency: `20`
- seed: `101`
- starting `Pi`: `0.99`

After-fix outputs:

- `benchmarks/reports/single_trait_test_data_pip_comparison_after_init_fix/run_summary.csv`
- `benchmarks/reports/single_trait_test_data_pip_comparison_after_init_fix/pi_trace_summary.csv`
- `benchmarks/reports/single_trait_test_data_pip_comparison_after_init_fix/annotated_class_summary.csv`

Summary:

| case | mean PIP before | mean PIP after |
| --- | ---: | ---: |
| full annotated BayesC short run | `0.06983` | `0.00669` |
| intercept-only annotated BayesC short run | `0.04992` | `0.01031` |

Trace summary after the fix:

| case | first saved inclusion | last saved inclusion | mean inclusion |
| --- | ---: | ---: | ---: |
| full annotated BayesC | `0.00631` | `0.00688` | `0.00679` |
| intercept-only annotated BayesC | `0.00692` | `0.01181` | `0.01041` |

The intercept-only short run now matches the expected sparse BayesC scale instead of spending the post-burnin samples near the inflated `0.05` baseline. The full annotated short run no longer has the high baseline either, but it now undershoots the earlier long full annotated run (`0.01946`), so full annotation-effect learning still needs enough chain length.
