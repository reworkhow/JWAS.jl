# Simulated Annotations Method Matrix Report

## Benchmark

Production-path benchmark on the packaged `simulated_annotations` dataset.

- dataset: `src/4.Datasets/data/simulated_annotations/`
- seeds: `100`, `110`
- `chain_length = 5000`
- `burnin = 1000`
- `output_samples_frequency = 10`
- starting `h2 = 0.5`

Output:

- `/tmp/simulated_annotations_method_matrix_100_110_len5000`
- `/tmp/simulated_annotations_method_matrix_100_110_len5000/cross_seed_summary.csv`

## Cross-Seed Results

| Variant | Marker corr | PIP corr | EBV corr | Annotation coeff corr | `pi` summary |
| --- | ---: | ---: | ---: | ---: | --- |
| `BayesC_dense` | `0.8502` | `0.6409` | `0.9813` | `NA` | scalar abs diff `0.3951` |
| `BayesC_fast_blocks_1` | `0.8512` | `0.6079` | `0.9845` | `NA` | scalar abs diff `0.3942` |
| `Annotated_BayesC_dense` | `0.9911` | `0.9932` | `0.9997` | `0.9956` | vector corr `0.99997` |
| `Annotated_BayesC_fast_blocks_1` | `0.3154` | `0.0092` | `0.8978` | `0.9534` | vector corr `0.0075` |
| `BayesR_dense` | `0.9760` | `0.8255` | `0.9995` | `NA` | vector corr `0.99992` |
| `BayesR_fast_blocks_1` | `0.9718` | `0.8032` | `0.9990` | `NA` | vector corr `0.99191` |
| `Annotated_BayesR_dense` | `0.9985` | `0.9982` | `0.9997` | `0.9844` | vector corr `0.999996` |
| `Annotated_BayesR_fast_blocks_1` | `0.9985` | `0.9982` | `0.9997` | `0.9844` | vector corr `0.999996` |

## Interpretation

On this packaged Jian-style simulated dataset:

- dense annotated BayesC is highly stable across seeds
- dense annotated BayesR is also highly stable across seeds
- plain BayesR is reasonably stable, but not as strong as the annotated version
- plain BayesC remains the weakest dense method on marker/PIP agreement
- `Annotated_BayesC_fast_blocks_1` is clearly unstable on this dataset
- `Annotated_BayesR_fast_blocks_1` is effectively identical to dense in this run

This report confirms that the packaged `simulated_annotations` dataset is a
useful benchmark fixture for comparing BayesC/BayesR annotation behavior across
seeds under a practical data regime.
