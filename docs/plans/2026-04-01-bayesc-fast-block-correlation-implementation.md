# BayesC Fast-Block Correlation Benchmark Implementation

## Goal

Run a same-seed production-path comparison for:

- dense `BayesC` vs `BayesC` with `fast_blocks=1`
- dense annotated `BayesC` vs annotated `BayesC` with `fast_blocks=1`

and quantify how close the realized posterior summaries remain.

## What Changed

Added one benchmark script:

- `benchmarks/bayesc_fast_block_correlation.jl`

Added two benchmark notes:

- `docs/plans/2026-04-01-bayesc-fast-block-correlation-design.md`
- `benchmarks/reports/2026-04-01-bayesc-fast-block-correlation-report.md`

## Benchmark Configuration

Final production-path run:

- data seed: `20260401`
- MCMC seeds: `2026,2027`
- `n_obs = 50`
- `n_markers = 100`
- `chain_length = 800`
- `burnin = 200`
- `output_samples_frequency = 100`
- synthetic annotated truth with the `stepwise_annotation_signal` conditional
  structure

Run command:

```bash
JWAS_BAYESC_CORR_CHAIN_LENGTH=800 \
JWAS_BAYESC_CORR_BURNIN=200 \
JWAS_BAYESC_CORR_OUTPUT_FREQ=100 \
JWAS_BAYESC_CORR_N_OBS=50 \
JWAS_BAYESC_CORR_N_MARKERS=100 \
JWAS_BAYESC_CORR_SEEDS=2026,2027 \
julia --project=. --startup-file=no \
  benchmarks/bayesc_fast_block_correlation.jl \
  /tmp/bayesc_fast_block_correlation_run
```

## Results Summary

The paired same-seed comparisons came back effectively identical for both
methods on this benchmark.

### Standard BayesC

- marker-effect correlation: `1.0`
- model-frequency correlation: `1.0`
- EBV correlation: `1.0`
- max marker-effect absolute difference for seed `2026`: about `4.5e-7`
- max EBV absolute difference for seed `2026`: about `4.49e-6`

### Annotated BayesC

- marker-effect correlation: `1.0`
- model-frequency correlation: `1.0`
- EBV correlation: `1.0`
- per-marker `pi_j` correlation: `1.0`
- annotation-coefficient correlation: `1.0`
- annotation-coefficient max absolute difference for seed `2026`: `0.0`

## Interpretation

On the current merged `master`, this benchmark does **not** reproduce a material
dense-vs-`fast_blocks=1` discrepancy for annotated BayesC. For these settings,
both the standard and annotated BayesC runs are effectively identical under
same-seed paired comparisons.

That means the earlier concern should be treated as unresolved for larger or
more extreme settings, not as an established current-master bug.
