# Remove Public Benchmark Page Implementation Record

## Summary

The fast-block benchmark page was removed from the public JWAS manual because
its target-scale times and speed ratio were extrapolated from smaller
calibration runs and could be mistaken for directly measured target-scale
performance guidance. The underlying engineering record remains available at
`benchmarks/reports/2026-02-26-fast-block-calibration-and-target-extrapolation.md`.

The archived report now distinguishes the evidence explicitly:

- wall-clock timings at `P = 100,000` and `P = 200,000` were measured;
- the estimates for `P = 2,000,000`, `chain_length = 2000`, and the approximately
  `213x` ratio were extrapolated for this calibration setup;
- the result is not a direct full-scale benchmark or a general performance
  guarantee; and
- the two Slurm calibration jobs ran on different CPU node types.

## Files Changed

- `docs/src/manual/benchmark.md` was moved out of the public manual to
  `benchmarks/reports/2026-02-26-fast-block-calibration-and-target-extrapolation.md`
  and reframed as an internal calibration and extrapolation report.
- `benchmarks/README.md` now describes the associated Slurm scripts as inputs to
  the archived calibration report instead of a public benchmark page.
- `docs/make.jl` no longer includes the `Benchmark` navigation entry.
- `docs/src/index.md` no longer includes the removed page in its source-page
  list.
- `docs/src/manual/block_bayesc.md` no longer links to or promotes the removed
  target-scale benchmark.
- `docs/src/manual/workflow.md` no longer links to or promotes the removed
  target-scale benchmark.
- `docs/plans/2026-07-11-remove-public-benchmark-page-plan.md` records the
  approved implementation plan.
- `docs/plans/2026-07-11-remove-public-benchmark-page-implementation.md` records
  the completed correction and verification evidence.

No production Julia source, public API, sampler behavior, tests, or benchmark
data changed. The measured calibration tables, extrapolation procedure, and
reported values were retained; only their framing and publication location
changed.

## Verification

The following commands and results correspond to the documentation correction
recorded at commit `1b6f9ce7`.

These checks were run from the repository worktree:

```bash
git diff --check 6b1cf84e..1b6f9ce7
test -f benchmarks/reports/2026-02-26-fast-block-calibration-and-target-extrapolation.md
test ! -e docs/src/manual/benchmark.md
if rg -n -i 'benchmark\.md|\[Benchmark\]|real cluster benchmark|target-scale benchmark' \
  docs/make.jl docs/src; then exit 1; fi
```

They completed successfully: the feature diff has no whitespace errors, the
internal report exists, the public source page is absent, and no stale public
references were found.

```bash
julia --project=docs --startup-file=no docs/make.jl
```

The documentation build exited successfully. Documenter's warning that it
could not auto-detect a deployment environment and skipped deployment is
expected for this local build. In the verification environment, `julia` was not
on `PATH`, so this equivalent command was run with Julia 1.11.7 via
`/Users/haocheng/.juliaup/bin/julia`.

```bash
test ! -e docs/build/manual/benchmark
if rg -n -i 'manual/benchmark|>Benchmark<|target-scale benchmark' docs/build; then
  exit 1
fi
```

The rendered site contains no `manual/benchmark` page, `Benchmark` navigation
label, or target-scale benchmark wording.

```bash
git status --short
git diff --stat 6b1cf84e..1b6f9ce7
git diff 6b1cf84e..1b6f9ce7 -- docs/make.jl docs/src benchmarks \
  docs/plans/2026-07-11-remove-public-benchmark-page-plan.md \
  docs/plans/2026-07-11-remove-public-benchmark-page-implementation.md
```

The cumulative feature diff is limited to the approved documentation
correction and its records. The local documentation environment and build
artifacts, `docs/Manifest.toml` and `docs/build/`, remain ignored and are not
tracked.
