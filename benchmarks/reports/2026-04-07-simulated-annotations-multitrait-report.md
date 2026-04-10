# Simulated Annotations Multi-Trait Benchmark Report

## Benchmark

Production-path benchmark on the packaged 2-trait extension of
`simulated_annotations`.

- benchmark script:
  [benchmarks/simulated_annotations_multitrait_comparison.jl](/Users/haocheng/Github/JWAS.jl/.worktrees/multitrait-annotated-bayesc/benchmarks/simulated_annotations_multitrait_comparison.jl)
- output directory:
  `/tmp/jwas_mt_sim_benchmark_20260407`
- seeds: `101`, `202`
- `chain_length = 2000`
- `burnin = 500`
- `output_samples_frequency = 20`
- benchmark warmup: enabled

The packaged 2-trait truth contains:

- `944` null markers (`00`)
- `6` trait-1-only markers (`10`)
- `6` trait-2-only markers (`01`)
- `8` shared markers (`11`)

So each trait has `14` active markers and there are `20` any-active markers.

## Method Summary

| Variant | Trait | `cor(y, EBV)` | Effect corr | Top-k recall | Any-active recall |
| --- | --- | ---: | ---: | ---: | ---: |
| `MT_BayesC` | `y1` | `0.7731` | `0.5057` | `0.2857` | `0.3000` |
| `MT_BayesC` | `y2` | `0.7353` | `0.7573` | `0.3571` | `0.3000` |
| `MT_Annotated_BayesC` | `y1` | `0.7607` | `0.4970` | `0.3571` | `0.4250` |
| `MT_Annotated_BayesC` | `y2` | `0.7146` | `0.7725` | `0.4643` | `0.4250` |
| `BayesC_y1` | `y1` | `0.7835` | `0.4723` | `0.2500` | `NA` |
| `Annotated_BayesC_y1` | `y1` | `0.7580` | `0.5403` | `0.3214` | `NA` |
| `BayesC_y2` | `y2` | `0.7279` | `0.7799` | `0.3571` | `NA` |
| `Annotated_BayesC_y2` | `y2` | `0.7069` | `0.8396` | `0.6071` | `NA` |
| `BayesR_y1` | `y1` | `0.8029` | `0.4752` | `0.2500` | `NA` |
| `BayesR_y2` | `y2` | `0.7468` | `0.8085` | `0.3929` | `NA` |

Single-trait family any-active recall, computed by combining the `y1` and `y2`
single-trait runs within each method family:

| Family | Any-active recall |
| --- | ---: |
| `BayesC_single` | `0.3000` |
| `Annotated_BayesC_single` | `0.6000` |
| `BayesR_single` | `0.3000` |

Joint shared-state posterior summary for the multi-trait runs, computed from the
saved marker-effect MCMC samples:

- for each marker, `P11 = Pr(δ = (1,1) | data)` is reconstructed as the
  frequency of simultaneous nonzero effects across the saved multi-trait MCMC
  samples
- a marker is declared shared when `P11 >= 0.5`
- single-trait runs are not included here because they do not define a joint
  `P11`

| Variant | Mean `P11` on true shared loci | Mean `P11` on nonshared loci | Declared shared (`P11 >= 0.5`) | True shared among declared | Shared precision | Shared recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `MT_BayesC` | `0.3700` | `0.0147` | `3.0` | `3.0` | `1.0000` | `0.3750` |
| `MT_Annotated_BayesC` | `0.4667` | `0.0118` | `5.5` | `3.0` | `0.5500` | `0.3750` |

## Interpretation

### 1. The new multi-trait annotated BayesC improves over conventional multi-trait BayesC on marker recovery.

Relative to `MT_BayesC`, `MT_Annotated_BayesC` improved:

- trait-1 top-k recall from `0.2857` to `0.3571`
- trait-2 top-k recall from `0.3571` to `0.4643`
- any-active recall from `0.3000` to `0.4250`
- trait-2 effect correlation from `0.7573` to `0.7725`

So on this fixture the annotation prior is helping the multi-trait BayesC model
prioritize active markers more effectively than the conventional multi-trait
baseline.

### 2. On the true joint `P11` metric, annotated multi-trait BayesC shifts more posterior mass toward the shared state but does not yet improve majority-threshold shared recall.

Relative to `MT_BayesC`, `MT_Annotated_BayesC` changed the shared-state
posterior in two useful ways:

- mean `P11` on the true shared loci increased from `0.3700` to `0.4667`
- mean `P11` on nonshared loci decreased from `0.0147` to `0.0118`

So the annotation model is pushing posterior mass in the intended direction.
But under the explicit declaration rule `P11 >= 0.5`, both methods recovered
the same number of true shared loci:

- `MT_BayesC`: `3.0` true shared declarations out of `8`
- `MT_Annotated_BayesC`: `3.0` true shared declarations out of `8`

The annotated model declared more shared markers overall (`5.5` vs `3.0`), so
its shared precision was lower at this threshold:

- `MT_BayesC`: precision `1.0000`
- `MT_Annotated_BayesC`: precision `0.5500`

So with a true joint posterior declaration rule, the current result is more
mixed than the earlier proxy ranking suggested.

### 3. The multi-trait annotated model does not yet dominate the strongest single-trait baselines on this scenario.

The strongest marker-recovery rows in this benchmark are the single-trait
annotated BayesC runs:

- `Annotated_BayesC_y1` had the best `y1` effect correlation: `0.5403`
- `Annotated_BayesC_y2` had the best `y2` effect correlation: `0.8396`
- `Annotated_BayesC_y2` had the best trait-specific top-k recall: `0.6071`
- the combined single-trait annotated BayesC family had the best any-active
  recall: `0.6000`

So the new multi-trait annotated BayesC is benchmark-positive relative to
ordinary multi-trait BayesC, but it is not yet the strongest marker-prioritizer
in this fixed simulated setting.

### 4. BayesR remained the strongest predictor on `cor(y, EBV)`.

Plain BayesR gave the highest prediction correlations on both traits:

- `BayesR_y1 = 0.8029`
- `BayesR_y2 = 0.7468`

Both multi-trait BayesC variants and both annotated BayesC single-trait runs
were lower on `cor(y, EBV)` in this benchmark.

### 5. The annotation signal was learned most clearly at step 1.

For multi-trait annotated BayesC, the step-1 `00` vs active model learned a
clear positive coefficient on `active_signal`:

- `step1_zero_vs_active`, `Annotation_1 = 1.0344`

The later steps were noisier across the two seeds:

- `step2_11_vs_singleton`, `Annotation_1 = 1.1073`
- `step3_10_vs_01`, `Annotation_2 = -2.5447`

So this benchmark shows strong evidence that the v1 implementation is using the
annotation information, but weaker evidence that steps 2 and 3 are already as
stable as step 1 under these settings.

## Conclusion

On the packaged 2-trait simulated annotation fixture:

1. Dense annotated 2-trait BayesC is working on the production path and beats
   conventional dense multi-trait BayesC on marker-prioritization metrics.
2. On the true joint shared-state posterior `P11`, it pushes more mass toward
   the shared state for true pleiotropic loci, but at the current `P11 >= 0.5`
   declaration rule it does not recover more true shared loci than
   conventional multi-trait BayesC.
3. It is not benchmark-dominant over the single-trait annotated BayesC
   baselines on this scenario.
4. Plain BayesR remains the best predictor by `cor(y, EBV)` in this report.

That is a useful benchmark result: the new method is not just runnable, it is
already competitive enough to improve the intended multi-trait baseline, while
still leaving clear room for future work on the step-2 and step-3 annotation
behavior.
