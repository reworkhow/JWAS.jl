# Multi-Trait Annotated BayesC Sampler Benchmark Report

## Benchmark

Production-path benchmark on the packaged 2-trait `simulated_annotations`
fixture using the updated benchmark harness in
`benchmarks/simulated_annotations_multitrait_comparison.jl`.

- output directory: `/tmp/jwas_mt_sampler_benchmark_20260410`
- seeds: `101`, `202`
- `chain_length = 2000`
- `burnin = 500`
- `output_samples_frequency = 20`
- benchmark warmup: disabled

The packaged truth contains:

- `944` null markers (`0`)
- `6` trait-1-only markers (`10`)
- `6` trait-2-only markers (encoded as `1` in `truth_mt.csv`)
- `8` shared markers (`11`)

So each trait has `14` active markers and there are `20` any-active markers.

## Method Summary

| Variant | Trait | Sampler | `cor(y, EBV)` | Effect corr | Top-k recall | Any-active recall | Runtime (s) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `MT_BayesC` | `y1` | `auto` | `0.7731` | `0.5057` | `0.2857` | `0.3000` | `12.12` |
| `MT_BayesC` | `y2` | `auto` | `0.7353` | `0.7573` | `0.3571` | `0.3000` | `12.12` |
| `MT_Annotated_BayesC_I` | `y1` | `:I` | `0.7607` | `0.4970` | `0.3571` | `0.4250` | `6.09` |
| `MT_Annotated_BayesC_I` | `y2` | `:I` | `0.7146` | `0.7725` | `0.4643` | `0.4250` | `6.09` |
| `MT_Annotated_BayesC_II` | `y1` | `:II` | `0.7550` | `0.5275` | `0.3929` | `0.5250` | `15.62` |
| `MT_Annotated_BayesC_II` | `y2` | `:II` | `0.7060` | `0.8378` | `0.6429` | `0.5250` | `15.62` |
| `BayesC_y1` | `y1` | `auto` | `0.7835` | `0.4723` | `0.2500` | `NA` | `1.14` |
| `Annotated_BayesC_y1` | `y1` | `auto` | `0.7580` | `0.5403` | `0.3214` | `NA` | `1.12` |
| `BayesC_y2` | `y2` | `auto` | `0.7279` | `0.7799` | `0.3571` | `NA` | `0.41` |
| `Annotated_BayesC_y2` | `y2` | `auto` | `0.7069` | `0.8396` | `0.6071` | `NA` | `0.75` |
| `BayesR_y1` | `y1` | `auto` | `0.8029` | `0.4752` | `0.2500` | `NA` | `1.02` |
| `BayesR_y2` | `y2` | `auto` | `0.7468` | `0.8085` | `0.3929` | `NA` | `0.90` |
| `Annotated_BayesR_y1` | `y1` | `auto` | `0.7805` | `0.5942` | `0.2500` | `NA` | `1.65` |
| `Annotated_BayesR_y2` | `y2` | `auto` | `0.7320` | `0.8620` | `0.5357` | `NA` | `1.49` |

Single-trait family any-active recall, computed by combining the `y1` and `y2`
single-trait runs within each family:

| Family | Any-active recall |
| --- | ---: |
| `Annotated_BayesC_single` | `0.6000` |
| `Annotated_BayesR_single` | `0.3750` |
| `BayesC_single` | `0.3000` |
| `BayesR_single` | `0.3000` |

## Pleiotropic Recovery

For the multi-trait methods, pleiotropic recovery uses the true joint posterior
shared-state score `P11 = Pr(delta = (1,1) | data)` reconstructed from the
saved multi-trait MCMC samples. Markers are ranked by `P11`, and the top `8`
markers are declared shared.

For the single-trait families, pleiotropy is a proxy based on overlap:

- detect top `14` trait-1 markers
- detect top `14` trait-2 markers
- declare shared markers from the overlap of those two detected sets

| Label | Shared score | Declared shared | True shared among declared | Shared precision | Shared recall | Shared F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `MT_BayesC` | `P11_topk` | `8.0` | `3.0` | `0.3750` | `0.3750` | `0.3750` |
| `MT_Annotated_BayesC_I` | `P11_topk` | `8.0` | `3.5` | `0.4375` | `0.4375` | `0.4375` |
| `MT_Annotated_BayesC_II` | `P11_topk` | `8.0` | `4.5` | `0.5625` | `0.5625` | `0.5625` |
| `Annotated_BayesC_single` | overlap top-k | `4.5` | `3.0` | `0.6071` | `0.3750` | `0.4333` |
| `Annotated_BayesR_single` | overlap top-k | `5.5` | `2.0` | `0.4583` | `0.2500` | `0.3068` |
| `BayesC_single` | overlap top-k | `2.5` | `2.5` | `1.0000` | `0.3125` | `0.4727` |
| `BayesR_single` | overlap top-k | `2.5` | `2.5` | `1.0000` | `0.3125` | `0.4727` |

Additional `P11` separation for the multi-trait methods:

| Label | Mean `P11` on true shared loci | Mean `P11` on nonshared loci |
| --- | ---: | ---: |
| `MT_BayesC` | `0.3700` | `0.0147` |
| `MT_Annotated_BayesC_I` | `0.4667` | `0.0118` |
| `MT_Annotated_BayesC_II` | `0.4925` | `0.0105` |

## Interpretation

### 1. Sampler II is the strongest pleiotropic detector among the tested methods.

`MT_Annotated_BayesC_II` is the best method in this report on the primary
pleiotropic-recovery metric:

- shared precision: `0.5625`
- shared recall: `0.5625`
- shared F1: `0.5625`

Relative to conventional multi-trait BayesC:

- true shared among declared top-8 improved from `3.0` to `4.5`
- shared recall improved from `0.3750` to `0.5625`
- mean `P11` on true shared loci increased from `0.3700` to `0.4925`
- mean `P11` on nonshared loci decreased from `0.0147` to `0.0105`

So the explicit sampler-II option is not just reachable; on this fixture it is
the best-performing multi-trait choice for the main pleiotropy target.

### 2. Sampler II also improves any-active and trait-specific recovery over sampler I, but it is slower.

Comparing the two annotated multi-trait samplers directly:

- any-active recall: `0.4250` with sampler I versus `0.5250` with sampler II
- `y1` top-k recall: `0.3571` with sampler I versus `0.3929` with sampler II
- `y2` top-k recall: `0.4643` with sampler I versus `0.6429` with sampler II

But the extra mixing freedom is not free:

- sampler I runtime mean: `6.09s`
- sampler II runtime mean: `15.62s`

So sampler II buys materially better marker recovery, especially for shared and
trait-2-active loci, at a noticeably higher runtime cost.

### 3. The genomic-prediction winners are still the single-trait BayesR baselines.

For `cor(y, EBV)`, plain BayesR remains best:

- `BayesR_y1 = 0.8029`
- `BayesR_y2 = 0.7468`

Both annotated multi-trait BayesC variants are below those values:

- sampler I: `0.7607` on `y1`, `0.7146` on `y2`
- sampler II: `0.7550` on `y1`, `0.7060` on `y2`

So the multi-trait annotated BayesC improvements here are about causal-marker
and pleiotropic recovery, not top-line genomic prediction accuracy.

### 4. Single-trait annotated BayesC remains a strong any-active baseline, but it does not match sampler II on shared recall.

`Annotated_BayesC_single` still has the best single-trait-family any-active
recall:

- `0.6000` versus `0.5250` for `MT_Annotated_BayesC_II`

But on shared recovery, the overlap proxy does not match the multi-trait
sampler-II joint posterior:

- `Annotated_BayesC_single` shared recall: `0.3750`
- `MT_Annotated_BayesC_II` shared recall: `0.5625`

That is the main benefit of the multi-trait approach in this comparison: it
recovers more true pleiotropic loci than the overlap of two single-trait runs.

### 5. The annotation signal is still clearest at step 1.

For both annotated multi-trait samplers, the step-1 active-vs-zero model learns
a strong positive coefficient on `Annotation_1`:

- sampler I: `1.0344`
- sampler II: `0.9419`

The later step-2 and step-3 coefficients are still more variable, especially
under sampler II. So the benchmark supports the feature strongly at the marker
recovery level, while still suggesting that later tree steps are noisier than
step 1 on this dataset.

## Conclusion

On the packaged `simulated_annotations` 2-trait fixture:

1. The new sampler-II option is useful. It is the best method in this report
   on the primary pleiotropic-recovery metric.
2. Multi-trait annotated BayesC with sampler II improves over both
   conventional multi-trait BayesC and annotated sampler I on shared recovery,
   any-active recovery, and trait-2 top-k recall.
3. The gain comes with a runtime cost and does not make annotated multi-trait
   BayesC the strongest genomic predictor by `cor(y, EBV)`.
4. Single-trait annotated BayesC remains a strong baseline for any-active
   detection, but the multi-trait sampler-II path now shows the clearest
   advantage on true pleiotropic recovery.
