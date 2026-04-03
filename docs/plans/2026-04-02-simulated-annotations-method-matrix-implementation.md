# Simulated Annotations Method Matrix Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

## Goal

Create and run a packaged benchmark script for cross-seed JWAS method
comparison on `simulated_annotations`.

## Tasks

### Task 1: Add the benchmark script

Create:

- `benchmarks/simulated_annotations_method_matrix.jl`

Requirements:

- use `JWAS.Datasets.dataset(..., dataset_name="simulated_annotations")`
- run the 8 benchmark variants
- accept optional environment overrides for seeds and MCMC settings
- write `cross_seed_summary.csv`

### Task 2: Run the benchmark

Run with the recommended defaults:

- seeds `100,110`
- `chain_length = 5000`
- `burnin = 1000`

Save output under `/tmp/`.

### Task 3: Save a report

Create:

- `benchmarks/reports/2026-04-02-simulated-annotations-method-matrix-report.md`

Include:

- dataset and MCMC settings
- output path
- cross-seed summary table
- brief interpretation

### Task 4: Report repository state

Call out:

- new benchmark script
- new report
- plan docs
- no production sampler changes
