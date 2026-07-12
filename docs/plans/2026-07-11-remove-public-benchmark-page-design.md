# Remove Public Benchmark Page Design

## Goal

Remove the current fast-block benchmark page from the public JWAS manual because
its target-scale headline is extrapolated from smaller calibration runs and can
be mistaken for a direct target-scale benchmark.

Preserve the measured calibration record under `benchmarks/reports/` for
engineering traceability.

## Rationale

The benchmark page contains real wall-clock measurements, but the reported
`N=50,000`, `P=2,000,000`, `chain_length=2000` times were not measured directly.
They were extrapolated from:

- marker counts of `100,000` and `200,000`;
- short chain or outer-iteration counts;
- two replicates;
- jobs that ran on different CPU node types.

That evidence is useful as an internal calibration record, but it is not strong
enough to support a prominent public target-scale performance comparison.

## Public Documentation Changes

- Remove `docs/src/manual/benchmark.md` from the public manual.
- Remove the `Benchmark` navigation entry from `docs/make.jl`.
- Remove the benchmark page entry from `docs/src/index.md`.
- Remove links that describe it as a real target-scale cluster benchmark from:
  - `docs/src/manual/block_bayesc.md`
  - `docs/src/manual/workflow.md`
- Keep the qualified theoretical complexity and unit-cost proxy discussion on
  the BayesR3 block-strategy page.
- Retain practical advice that users should benchmark their own workloads and
  block choices.

## Internal Record

Move the benchmark content to:

`benchmarks/reports/2026-02-26-fast-block-calibration-and-target-extrapolation.md`

The internal report will state prominently that:

- the calibration timings were measured;
- the target-scale times and speed ratio were extrapolated;
- the result is not a direct full-scale benchmark or a general performance
  guarantee;
- the two benchmark jobs ran on different CPU node types.

## Scope

This is a documentation-only correction. It changes no Julia source code,
sampler behavior, public API, benchmark data, or test behavior.

## Verification

1. Search `docs/make.jl` and `docs/src` for stale links to `benchmark.md`.
2. Confirm the internal report exists under `benchmarks/reports/`.
3. Run `git diff --check`.
4. Run `/Users/haocheng/.juliaup/bin/julia --project=docs --startup-file=no docs/make.jl`.
5. Confirm the rendered navigation no longer contains the removed Benchmark
   page.
6. Add a matching implementation record under `docs/plans/`.
