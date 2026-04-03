# Simulated Annotations Method Matrix Design

## Goal

Add a reproducible benchmark-style regression test that runs the main dense and
`fast_blocks=1` BayesC/BayesR variants on the packaged
`simulated_annotations` dataset and summarizes cross-seed stability.

## Scope

Benchmark cases:

- `BayesC_dense`
- `BayesC_fast_blocks_1`
- `Annotated_BayesC_dense`
- `Annotated_BayesC_fast_blocks_1`
- `BayesR_dense`
- `BayesR_fast_blocks_1`
- `Annotated_BayesR_dense`
- `Annotated_BayesR_fast_blocks_1`

The benchmark should:

- resolve the packaged dataset through `JWAS.Datasets`
- run two seeds on the same saved dataset
- summarize marker, PIP, EBV, annotation-coefficient, and `pi` agreement
- write a CSV summary to an output directory

## Why A Benchmark Script Instead Of A Unit Test

This check requires 16 production MCMC runs. That is too heavy for the regular
unit suite and is better treated as a reproducible benchmark/regression script.

## Recommended Defaults

- dataset: `simulated_annotations`
- seeds: `100,110`
- `chain_length = 5000`
- `burnin = 1000`
- `output_samples_frequency = 10`
- starting `h2 = 0.5`

## Output

- benchmark script in `benchmarks/`
- one markdown report in `benchmarks/reports/`
