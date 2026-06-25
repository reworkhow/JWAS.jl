# Multi-Trait Annotated BayesR Simulated Annotations Smoke Test

Date: 2026-06-24

This report records a focused smoke simulation for the new multi-trait annotated
BayesR path using the bundled `simulated_annotations` genotype and annotation
data.

## Command

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=mt_annotated_bayesr \
JWAS_SIMULATED_MT_SEEDS=101 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=300 \
JWAS_SIMULATED_MT_BURNIN=100 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=20 \
JWAS_SIMULATED_MT_WARMUP=false \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  benchmarks/out/mt_annotated_bayesr_smoke_20260624
```

## Data

- Genotypes: `src/4.Datasets/data/simulated_annotations/genotypes.csv`
- Phenotypes: `src/4.Datasets/data/simulated_annotations/phenotypes_mt.csv`
- Annotations: `src/4.Datasets/data/simulated_annotations/annotations_mt.csv`
- Truth: `src/4.Datasets/data/simulated_annotations/truth_mt.csv`
- Marker count: 964
- Individual count: 400

## Checks

The benchmark harness normalizes marker IDs, joins marker outputs to truth by
`marker_id`, and asserts:

- output and truth have the same marker count
- marker key sets match
- no duplicate marker keys are present
- the joined row count matches the truth row count

The smoke run produced 964 joined marker rows for each trait and 964 rows in
`posterior_joint_state_probabilities.csv`.

## Results

Only one seed and 10 saved posterior samples were used, so these numbers test
that the production path runs and summarizes correctly. They should not be read
as convergence or method-quality estimates.

| Variant | Trait | Runtime (s) | EBV correlation | Marker effect correlation | PIP gap | Top-k recall |
|---|---:|---:|---:|---:|---:|---:|
| MT_Annotated_BayesR | y1 | 13.37 | 0.798 | 0.378 | 0.301 | 0.071 |
| MT_Annotated_BayesR | y2 | 13.37 | 0.758 | 0.646 | 0.296 | 0.143 |

Shared-effect posterior summary:

| Shared truth count | Declared shared count | True shared declared | Precision | Recall | F1 |
|---:|---:|---:|---:|---:|---:|
| 8 | 8 | 1 | 0.125 | 0.125 | 0.125 |

The annotation output included the expected seven multi-trait BayesR probit
steps:

- `zero_vs_active`
- `11_vs_singleton`
- `10_vs_01`
- `trait1_medium_or_large_vs_small`
- `trait1_large_vs_medium`
- `trait2_medium_or_large_vs_small`
- `trait2_large_vs_medium`

## Interpretation

The focused simulation verifies that annotated multi-trait BayesR can run on the
bundled simulated-annotations data, write marker-effect samples for both traits,
estimate annotation coefficients for the seven nested steps, and pass the
existing marker-key-aligned summary checks.

The short chain is intentionally a smoke test. A production comparison should
use longer chains, multiple seeds, and BayesC or single-trait BayesR baselines.
