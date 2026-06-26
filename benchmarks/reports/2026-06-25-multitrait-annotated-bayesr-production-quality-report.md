# Multi-Trait Annotated BayesR Production Quality Benchmark

## Goal

Evaluate whether the current multi-trait annotated BayesR implementation is ready for causal-variant discovery, not only prediction. The benchmark uses the production JWAS analysis path on `data/simulated_annotations`, joins marker outputs to truth by explicit `marker_id`, and summarizes causal discovery accuracy and cross-seed stability.

## Commands

Baseline comparison across the new method and BayesC/BayesR baselines:

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=bayesr_quality \
JWAS_SIMULATED_MT_SEEDS=101,202,303,404,505 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=10000 \
JWAS_SIMULATED_MT_BURNIN=2000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=20 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  benchmarks/out/mt_annotated_bayesr_quality_20260625
```

Long-chain confirmatory run for only `MT_Annotated_BayesR`:

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=mt_annotated_bayesr \
JWAS_SIMULATED_MT_SEEDS=101,202,303,404,505 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=50000 \
JWAS_SIMULATED_MT_BURNIN=10000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=50 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  benchmarks/out/mt_annotated_bayesr_quality_long_20260625
```

Marker alignment was checked after both runs: 50 baseline joined marker tables and 15 long-chain joined marker tables each had 964 unique canonical marker IDs with matching key sets.

## Causal Discovery

Top-k uses k equal to the number of true active markers for the trait. There are 14 active markers per trait, 20 any-active markers, and 8 shared active markers.

| Run | Method | Trait | AP | Top-k precision | PIP gap | Mean active PIP | Mean inactive PIP | EBV corr |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 10k | MT_Annotated_BayesC_I | y1 | 0.503 | 0.457 | 0.400 | 0.423 | 0.023 | 0.755 |
| 10k | MT_Annotated_BayesC_I | y2 | 0.580 | 0.486 | 0.453 | 0.471 | 0.019 | 0.709 |
| 10k | MT_Annotated_BayesR | y1 | 0.306 | 0.271 | 0.449 | 0.762 | 0.313 | 0.782 |
| 10k | MT_Annotated_BayesR | y2 | 0.470 | 0.429 | 0.490 | 0.683 | 0.193 | 0.717 |
| 10k | BayesR_y1 | y1 | 0.246 | 0.243 | 0.125 | 0.446 | 0.321 | 0.805 |
| 10k | BayesR_y2 | y2 | 0.412 | 0.414 | 0.284 | 0.430 | 0.146 | 0.752 |
| 10k | Annotated_BayesR_y1 | y1 | 0.295 | 0.286 | 0.496 | 0.760 | 0.264 | 0.794 |
| 10k | Annotated_BayesR_y2 | y2 | 0.625 | 0.571 | 0.541 | 0.592 | 0.051 | 0.723 |
| 50k | MT_Annotated_BayesR | y1 | 0.464 | 0.429 | 0.499 | 0.733 | 0.234 | 0.785 |
| 50k | MT_Annotated_BayesR | y2 | 0.598 | 0.543 | 0.572 | 0.626 | 0.055 | 0.720 |

The 10k run is not enough for the new method. It predicts reasonably, but its causal discovery is unstable and assigns too much PIP to inactive markers. The 50k run improves substantially: trait 2 becomes competitive with BayesC and annotated single-trait BayesR, while trait 1 improves but still carries high inactive-marker PIP.

## Stability

| Run | Method | Trait | PIP corr | Top-k Jaccard | Effect corr |
| --- | --- | --- | ---: | ---: | ---: |
| 10k | MT_Annotated_BayesC_I | y1 | 0.642 | 0.272 | 0.795 |
| 10k | MT_Annotated_BayesC_I | y2 | 0.711 | 0.391 | 0.909 |
| 10k | MT_Annotated_BayesC_I | any_active | 0.689 | 0.359 | 0.852 |
| 10k | MT_Annotated_BayesR | y1 | 0.276 | 0.124 | 0.833 |
| 10k | MT_Annotated_BayesR | y2 | 0.377 | 0.232 | 0.915 |
| 10k | MT_Annotated_BayesR | any_active | 0.281 | 0.146 | 0.860 |
| 50k | MT_Annotated_BayesR | y1 | 0.622 | 0.366 | 0.886 |
| 50k | MT_Annotated_BayesR | y2 | 0.685 | 0.530 | 0.947 |
| 50k | MT_Annotated_BayesR | any_active | 0.623 | 0.382 | 0.914 |

Longer chains largely resolve the most obvious stability problem. The 50k MT annotated BayesR PIP correlations are close to the 10k BayesC baseline, and top-k overlap is better than BayesC for y2 and any-active.

## Pleiotropy

| Run | Method | Shared precision | Shared recall | Shared F1 | True shared P11 | Non-shared P11 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 10k | MT_Annotated_BayesC_I | 0.575 | 0.575 | 0.575 | 0.558 | 0.013 |
| 10k | MT_Annotated_BayesR | 0.500 | 0.500 | 0.500 | 0.823 | 0.173 |
| 50k | MT_Annotated_BayesR | 0.550 | 0.550 | 0.550 | 0.781 | 0.055 |

The long-chain BayesR run improves shared-marker discrimination and greatly reduces non-shared P11, but BayesC still has cleaner non-shared P11 in this simulation.

## Recommendation

The implementation is technically usable and the method is promising, but I would not call the current default 10k-scale run production-ready for causal-variant discovery. The new method needs longer chains for stable PIPs. For this dataset, 50k iterations with 10k burn-in gives much more defensible results, but trait 1 still shows inflated inactive-marker PIP.

Before using this as a production causal-discovery method, the next checks should be:

- run a longer multi-method comparison at 50k+ for BayesC and all BayesR baselines, not only the new method
- test the proposed active-only `G` update variant against the current all-marker latent `beta` update
- inspect calibration of BayesR magnitude probabilities and PIP thresholds, especially for trait 1
- repeat on additional simulated architectures with different shared/singleton ratios and annotation strengths

Raw outputs are under `benchmarks/out/mt_annotated_bayesr_quality_20260625` and `benchmarks/out/mt_annotated_bayesr_quality_long_20260625`.
