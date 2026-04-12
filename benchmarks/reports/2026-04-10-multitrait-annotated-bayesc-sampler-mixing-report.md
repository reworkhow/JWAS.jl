# Multi-Trait Annotated BayesC Sampler Mixing Follow-Up Report

## Benchmark

Focused production-path rerun on the packaged 2-trait `simulated_annotations`
fixture using the benchmark harness in
`benchmarks/simulated_annotations_multitrait_comparison.jl`.

This follow-up keeps only the three multi-trait methods needed for the sampler
question:

- `MT_BayesC`
- `MT_Annotated_BayesC_I`
- `MT_Annotated_BayesC_II`

The run used:

- output directory: `/tmp/jwas_mt_sampler_mixing_20260410`
- seeds: `101`, `202`, `303`, `404`
- `chain_length = 4000`
- `burnin = 1000`
- `output_samples_frequency = 40`
- benchmark warmup: disabled

The packaged truth contains:

- `944` null markers
- `6` trait-1-only markers
- `6` trait-2-only markers
- `8` shared markers

So each trait has `14` active markers and there are `20` any-active markers.

## Exact-Posterior Backstop

Before interpreting the production rerun, the sampler comparison now has a
tiny-case regression in `test/unit/test_multitrait_mcmc.jl`.

That test constructs a one-marker, two-trait annotated BayesC case, computes
the exact posterior over the four joint states `00`, `10`, `01`, and `11`,
then checks that:

- sampler I empirical state frequencies match the exact posterior within Monte
  Carlo tolerance
- sampler II empirical state frequencies match the same exact posterior within
  Monte Carlo tolerance

So the remaining production benchmark question is finite-chain mixing and
efficiency, not whether the two samplers target different posteriors.

## Method Summary

| Variant | Trait | Sampler | `cor(y, EBV)` mean | Effect corr mean | Top-k recall mean | Any-active recall mean | Runtime mean (s) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `MT_BayesC` | `y1` | `auto` | `0.7741` | `0.5055` | `0.3036` | `0.3625` | `19.19` |
| `MT_BayesC` | `y2` | `auto` | `0.7294` | `0.7804` | `0.4464` | `0.3625` | `19.19` |
| `MT_Annotated_BayesC_I` | `y1` | `:I` | `0.7546` | `0.5341` | `0.3929` | `0.4875` | `12.99` |
| `MT_Annotated_BayesC_I` | `y2` | `:I` | `0.7111` | `0.8029` | `0.5179` | `0.4875` | `12.99` |
| `MT_Annotated_BayesC_II` | `y1` | `:II` | `0.7611` | `0.5727` | `0.4286` | `0.5000` | `35.18` |
| `MT_Annotated_BayesC_II` | `y2` | `:II` | `0.7106` | `0.8122` | `0.5357` | `0.5000` | `35.18` |

Selected standard deviations across seeds:

| Variant | Runtime sd (s) | `y1` top-k sd | `y2` top-k sd | Any-active recall sd |
| --- | ---: | ---: | ---: | ---: |
| `MT_BayesC` | `15.45` | `0.0357` | `0.0357` | `0.0250` |
| `MT_Annotated_BayesC_I` | `4.39` | `0.0412` | `0.1071` | `0.1109` |
| `MT_Annotated_BayesC_II` | `1.97` | `0.1010` | `0.0714` | `0.0408` |

## Pleiotropic Recovery

Pleiotropic recovery uses the true joint posterior shared-state score:

- `P11 = Pr(delta = (1,1) | data)`

For each seed and method:

- reconstruct `P11` from the saved multi-trait MCMC samples
- rank markers by `P11`
- declare the top `8` markers as shared

### Summary across seeds

| Variant | True shared among top 8 mean | Shared precision mean | Shared recall mean | Shared F1 mean | Mean `P11` on true shared | Mean `P11` on nonshared |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `MT_BayesC` | `3.00` | `0.3750` | `0.3750` | `0.3750` | `0.3696` | `0.0149` |
| `MT_Annotated_BayesC_I` | `3.00` | `0.3750` | `0.3750` | `0.3750` | `0.4867` | `0.0131` |
| `MT_Annotated_BayesC_II` | `4.25` | `0.5313` | `0.5313` | `0.5313` | `0.5638` | `0.0147` |

### Variability across seeds

| Variant | Shared precision sd | Shared recall sd | Shared F1 sd | `P11` true-shared sd | `P11` nonshared sd |
| --- | ---: | ---: | ---: | ---: | ---: |
| `MT_BayesC` | `0.0000` | `0.0000` | `0.0000` | `0.0071` | `0.0044` |
| `MT_Annotated_BayesC_I` | `0.0000` | `0.0000` | `0.0000` | `0.0713` | `0.0066` |
| `MT_Annotated_BayesC_II` | `0.1573` | `0.1573` | `0.1573` | `0.1115` | `0.0037` |

Per-seed true shared counts inside the top-8 shared set:

| Variant | Seed 101 | Seed 202 | Seed 303 | Seed 404 |
| --- | ---: | ---: | ---: | ---: |
| `MT_BayesC` | `3` | `3` | `3` | `3` |
| `MT_Annotated_BayesC_I` | `3` | `3` | `3` | `3` |
| `MT_Annotated_BayesC_II` | `4` | `4` | `6` | `3` |

## Interpretation

### 1. Sampler I and sampler II do not look different because of different target posteriors.

The exact tiny-case regression supports a common target posterior for both
samplers. So the production differences below should be read as finite-chain
mixing and efficiency differences.

### 2. After longer chains, sampler I no longer shows a pleiotropic advantage over conventional multi-trait BayesC.

On the primary shared-recovery metric:

- `MT_BayesC` shared recall = `0.3750`
- `MT_Annotated_BayesC_I` shared recall = `0.3750`

Both methods recover `3` true shared markers on average inside the top-8
shared set.

Sampler I still raises mean `P11` on true shared loci:

- `0.3696` for `MT_BayesC`
- `0.4867` for `MT_Annotated_BayesC_I`

but that increase did not translate into better top-8 shared recovery at this
longer chain setting.

### 3. Sampler II still has the strongest mean pleiotropic recovery, but it remains variable across seeds.

Sampler II is still best on the primary metric:

- shared precision/recall/F1 mean = `0.5313`
- true shared among declared top-8 mean = `4.25`
- mean `P11` on true shared loci = `0.5638`

This is an improvement over both:

- `0.3750` shared recall for `MT_BayesC`
- `0.3750` shared recall for `MT_Annotated_BayesC_I`

But the advantage is not fully stable across seeds:

- seed 303 recovered `6` shared markers
- seed 404 recovered only `3`
- shared recall sd = `0.1573`

So sampler II still looks stronger on average, but the benchmark does not yet
support calling that advantage fully stable.

### 4. Sampler II pays a large runtime cost for its shared-state advantage.

Mean runtimes:

- `MT_Annotated_BayesC_I`: `12.99s`
- `MT_BayesC`: `19.19s`
- `MT_Annotated_BayesC_II`: `35.18s`

So sampler II costs about:

- `2.71x` sampler I runtime
- `1.83x` conventional multi-trait BayesC runtime

That cost is material. Any recommendation to prefer sampler II should be tied
to stronger shared-state recovery being the main objective.

### 5. Secondary recovery metrics still favor the annotated methods over conventional multi-trait BayesC.

Any-active recall:

- `MT_BayesC`: `0.3625`
- `MT_Annotated_BayesC_I`: `0.4875`
- `MT_Annotated_BayesC_II`: `0.5000`

Per-trait top-k recall also remains better on average for the annotated
methods, especially on `y2`.

So even though sampler I does not improve top-8 shared recovery over
conventional multi-trait BayesC, the annotation model is still helping on
broader active-marker recovery.

### 6. Genomic prediction remains a secondary strength, not the main win.

Prediction accuracy for genomic prediction, measured as `cor(y, EBV)`, does
not clearly favor the annotated samplers:

- `MT_BayesC`: `0.7741` on `y1`, `0.7294` on `y2`
- `MT_Annotated_BayesC_I`: `0.7546` on `y1`, `0.7111` on `y2`
- `MT_Annotated_BayesC_II`: `0.7611` on `y1`, `0.7106` on `y2`

So the practical case for annotated multi-trait BayesC remains marker recovery
and pleiotropy, not top-line prediction accuracy.

## Conclusion

On the packaged `simulated_annotations` 2-trait fixture, after a longer-chain
focused rerun:

1. sampler I and sampler II should still be treated as the same-target
   samplers, because the exact tiny-case regression supports that conclusion
2. sampler I does not show a shared-recovery gain over conventional
   multi-trait BayesC at this longer chain setting
3. sampler II still gives the best mean pleiotropic recovery, but with
   substantial runtime cost and nontrivial between-seed variability
4. the remaining sampler question is practical mixing efficiency, not
   posterior correctness

The next careful step, if needed, is not another short benchmark. It is either:

- even longer chains with additional seeds, or
- explicit trace and ESS-style diagnostics on `P11`-driven shared-state
  summaries

before making a stronger production recommendation between sampler I and
sampler II.
