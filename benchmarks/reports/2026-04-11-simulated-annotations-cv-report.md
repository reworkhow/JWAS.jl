# Simulated Annotations Cross-Validation Report

## Benchmark Setup

Production K-fold cross-validation on the packaged `simulated_annotations`
fixture using the production harness in
`benchmarks/simulated_annotations_multitrait_comparison.jl`.

Configuration:

- benchmark mode: `cv`
- folds: `5`
- seeds: `101`, `202`
- chain length: `1500`
- burnin: `500`
- output samples frequency: `50`
- warmup: enabled
- output directory: `/tmp/jwas_simulated_annotations_cv_full_20260411`

Cross-validation rule:

- the same individual-level folds were used for every method within a seed
- for multi-trait methods, both `y1` and `y2` were hidden in the held-out fold
- for single-trait methods, only the target trait was hidden in the held-out
  fold
- genomic prediction metrics were computed only on held-out individuals

Primary CV metric:

- held-out prediction correlation `cor(y, EBV)`

Secondary CV metric:

- held-out RMSE

## Method Matrix

Multi-trait methods:

- `MT_BayesC`
- `MT_BayesC_I`
- `MT_BayesC_II`
- `MT_Annotated_BayesC_I`
- `MT_Annotated_BayesC_II`
- `MT_EmptyAnnotated_BayesC_I`
- `MT_EmptyAnnotated_BayesC_II`

Single-trait methods:

- `BayesC_single`
- `Annotated_BayesC_single`
- `BayesR_single`
- `Annotated_BayesR_single`

`MT_EmptyAnnotated_BayesC_*` remains a diagnostic control, not a recommended
user-facing analysis. It is included here because the requested multi-trait
comparison covered all current multi-trait models.

## Multi-Trait Family Summary

Trait-mean held-out prediction summary across `y1` and `y2`:

| Family | Mean held-out `cor(y, EBV)` | Mean held-out RMSE | Mean runtime per fold (s) |
| --- | ---: | ---: | ---: |
| `MT_Annotated_BayesC_II` | `0.6522` | `4.4700` | `9.40` |
| `MT_EmptyAnnotated_BayesC_II` | `0.6446` | `4.3642` | `9.01` |
| `MT_BayesC_II` | `0.6423` | `4.6572` | `8.80` |
| `MT_BayesC` | `0.6397` | `5.0015` | `2.22` |
| `MT_BayesC_I` | `0.6397` | `5.0015` | `2.26` |
| `MT_EmptyAnnotated_BayesC_I` | `0.6348` | `5.0477` | `2.30` |
| `MT_Annotated_BayesC_I` | `0.6338` | `4.5087` | `2.55` |

Main reading:

1. `MT_BayesC` and `MT_BayesC_I` are identical on held-out prediction here.
   That is the expected result because `:auto` dispatch uses sampler I on this
   plain multi-trait setup.
2. For both the plain and real-annotated multi-trait families, sampler II
   improves mean held-out correlation relative to sampler I, but the runtime
   cost is about `3.7x` to `4.0x`.
3. `MT_Annotated_BayesC_II` is the best multi-trait family on the primary CV
   correlation metric.

## Single-Trait Family Summary

Trait-mean held-out prediction summary across `y1` and `y2`:

| Family | Mean held-out `cor(y, EBV)` | Mean held-out RMSE | Mean runtime per fold (s) |
| --- | ---: | ---: | ---: |
| `BayesR_single` | `0.6497` | `4.4392` | `0.44` |
| `Annotated_BayesR_single` | `0.6484` | `4.3353` | `0.84` |
| `Annotated_BayesC_single` | `0.6469` | `4.6178` | `0.41` |
| `BayesC_single` | `0.6424` | `4.5294` | `0.19` |

Main reading:

1. Single-trait `BayesR` remains the strongest single-trait family on held-out
   prediction correlation.
2. Annotation improves the single-trait BayesC family mean
   (`0.6469` vs `0.6424`), but does not improve the single-trait BayesR family
   mean (`0.6484` vs `0.6497`).
3. Single-trait methods are much faster than the multi-trait sampler-II
   variants.

## Trait-Level Leaders

### `y1`

Top held-out correlation:

| Variant | Held-out `cor(y1, EBV_y1)` | Held-out RMSE | Mean runtime per fold (s) |
| --- | ---: | ---: | ---: |
| `BayesR_y1` | `0.6875` | `4.7027` | `0.44` |
| `MT_Annotated_BayesC_II` | `0.6875` | `4.9675` | `9.40` |
| `MT_EmptyAnnotated_BayesC_II` | `0.6865` | `4.0448` | `9.01` |
| `MT_EmptyAnnotated_BayesC_I` | `0.6819` | `4.6828` | `2.30` |
| `MT_BayesC_II` | `0.6801` | `4.4579` | `8.80` |

`BayesR_y1` and `MT_Annotated_BayesC_II` are effectively tied on held-out
correlation.

### `y2`

Top held-out correlation:

| Variant | Held-out `cor(y2, EBV_y2)` | Held-out RMSE | Mean runtime per fold (s) |
| --- | ---: | ---: | ---: |
| `Annotated_BayesC_y2` | `0.6232` | `4.4506` | `0.41` |
| `Annotated_BayesR_y2` | `0.6198` | `3.7568` | `0.83` |
| `MT_Annotated_BayesC_II` | `0.6168` | `3.9726` | `9.40` |
| `BayesR_y2` | `0.6118` | `4.1758` | `0.44` |
| `BayesC_y2` | `0.6076` | `4.6934` | `0.19` |

For `y2`, the strongest held-out correlation came from the single-trait
annotated BayesC run rather than a multi-trait method.

## Main Conclusions

1. Cross-validation changes the story relative to the in-sample EBV benchmark.
   The strongest family-level held-out predictor in this run was
   `MT_Annotated_BayesC_II`, not plain multi-trait BayesC.
2. The gain from sampler II is real on out-of-fold prediction, but it is
   expensive. For the real annotated multi-trait family:
   - sampler I: `0.6338` trait-mean held-out correlation, `2.55s` per fold
   - sampler II: `0.6522` trait-mean held-out correlation, `9.40s` per fold
3. Plain multi-trait BayesC gains only modestly from sampler II on this CV
   benchmark:
   - `MT_BayesC` / `MT_BayesC_I`: `0.6397`
   - `MT_BayesC_II`: `0.6423`
4. Single-trait BayesR remains highly competitive for prediction. It is the
   best single-trait family on mean held-out correlation and essentially tied
   for the best `y1` trait-level result.
5. The empty-annotated controls are prediction-competitive, especially with
   sampler II. That is interesting diagnostically, but it should not be
   over-interpreted as a preferred scientific model without further study.
