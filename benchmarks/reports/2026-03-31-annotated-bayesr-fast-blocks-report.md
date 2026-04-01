# Annotated BayesR Fast-Blocks Report

## Goal

Enable `fast_blocks` for dense single-trait annotated BayesR on the production
JWAS path and benchmark how the block sampler behaves relative to dense
annotated BayesR.

The implementation question was narrow:

- does the existing BayesR block kernel already consume annotation-driven
  marker-specific priors correctly?

The scientific question was separate:

- once enabled, how close are block annotated BayesR runs to dense annotated
  BayesR on a production benchmark?

## Implementation Scope

The block kernel itself was already annotation-ready. The enabling changes were:

- remove the validation guard that rejected annotated BayesR with
  `fast_blocks != false`
- replace the unit test that expected failure with a real annotated block run
- extend the annotated BayesR production benchmark with:
  - `Annotated_BayesR_fast_blocks_default`
  - `Annotated_BayesR_fast_blocks_1`
- document dense block support in the annotated BayesR manual

No BayesR block-kernel formulas were changed in this slice.

## Protocol

Production benchmark script:

- `benchmarks/annotated_bayesr_comparison.jl`

Output directory:

- `/tmp/annotated_bayesr_fast_blocks_benchmark_20260331`

Scenario:

- `stepwise_annotation_signal`

Common settings:

- single trait
- dense storage
- `n_obs = 200`
- `n_markers = 1000`
- target heritability `0.45`
- `chain_length = 10000`
- `burnin = 2000`
- `output_samples_frequency = 10`
- seeds `2026, 2027, 2028, 2029, 2030`

Method variants:

- `BayesR`
- `Annotated_BayesR`
- `Annotated_BayesR_fast_blocks_1`
- `Annotated_BayesR_fast_blocks_default`
- `Annotated_BayesC`

The stepwise-signal truth was the same as the earlier annotated BayesR report:

- baseline:
  - `P(delta > 1) = 0.10`
  - `P(delta > 2 | delta > 1) = 0.20`
  - `P(delta > 3 | delta > 2) = 0.20`
- enriched:
  - `P(delta > 1) = 0.30`
  - `P(delta > 2 | delta > 1) = 0.60`
  - `P(delta > 3 | delta > 2) = 0.60`

## Results

Mean results across five seeds:

| Method | fast_blocks | block size | outer chain | cor(y, EBV) | mean PIP causal | mean PIP null | enriched PIP | baseline PIP | top-causal recall | seconds |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `Annotated_BayesR` | `false` | `1` | `10000` | `0.7409` | `0.1027` | `0.0606` | `0.1388` | `0.0426` | `0.2986` | `6.69` |
| `Annotated_BayesR_fast_blocks_1` | `1` | `1` | `10000` | `0.7871` | `0.1156` | `0.0672` | `0.1632` | `0.0445` | `0.2944` | `7.70` |
| `Annotated_BayesR_fast_blocks_default` | `default` | `14` | `714` | `0.7792` | `0.1043` | `0.0556` | `0.1548` | `0.0319` | `0.2583` | `2.63` |

Reference rows from the same benchmark:

| Method | cor(y, EBV) | top-causal recall | seconds |
| --- | ---: | ---: | ---: |
| `BayesR` | `0.9164` | `0.2222` | `6.47` |
| `Annotated_BayesC` | `0.8595` | `0.2472` | `4.08` |

## Dense vs Block Deltas

Mean absolute differences versus dense annotated BayesR across matched seeds:

### `Annotated_BayesR_fast_blocks_1`

- `cor(y, EBV)`: `0.0520`
- residual variance: `0.0788`
- mean PIP causal: `0.0495`
- mean PIP null: `0.0379`
- enriched mean PIP: `0.0786`
- top-causal recall: `0.0097`
- `sigmaSq`: `0.2509`

### `Annotated_BayesR_fast_blocks_default`

- `cor(y, EBV)`: `0.0524`
- residual variance: `0.0938`
- mean PIP causal: `0.0732`
- mean PIP null: `0.0488`
- enriched mean PIP: `0.1294`
- top-causal recall: `0.0403`
- `sigmaSq`: `0.1283`

## Annotation-Coefficient Behavior

Dense annotated BayesR step means:

- step 1 `Annotation_1`: `+1.3349`
- step 2 `Annotation_1`: `+0.2365`
- step 3 `Annotation_1`: `+0.2518`

`fast_blocks=1` step means:

- step 1 `Annotation_1`: `+1.0823`
- step 2 `Annotation_1`: `+0.3156`
- step 3 `Annotation_1`: `+0.4830`

Default fast blocks step means:

- step 1 `Annotation_1`: `+0.9192`
- step 2 `Annotation_1`: `+0.0696`
- step 3 `Annotation_1`: `+1.1223`

So the block sampler still learns the informative annotation direction, but the
deeper step coefficients move materially relative to dense annotated BayesR.

## Interpretation

The engineering conclusion is straightforward:

- annotated BayesR block runs are now supported and execute correctly on the
  production path

The scientific conclusion is more cautious:

- block annotated BayesR is not a dense-parity mode
- even `fast_blocks=1` does not reproduce dense annotated BayesR closely on
  this benchmark
- default `fast_blocks=true` is much faster than dense annotated BayesR
  (`2.63s` vs `6.69s` mean runtime), but it shifts posterior summaries
  materially

This differs from the earlier ordinary BayesR fast-block work, where
`fast_blocks=1` was effectively a dense-parity setting. Once the annotation
layer is active, the dense/block differences are amplified enough to show up in
PIP, `sigmaSq`, and EBV correlation summaries.

So the correct user-facing description is:

- annotated BayesR with `fast_blocks` is supported
- it should be treated as an approximation to dense annotated BayesR, not as an
  equivalent mode

## Files

- benchmark script:
  - `benchmarks/annotated_bayesr_comparison.jl`
- raw outputs:
  - `/tmp/annotated_bayesr_fast_blocks_benchmark_20260331`
