# Independent-Blocks Fast-Block Parallelism Report

**Date:** 2026-04-15

**Feature branch:** `feature/independent-blocks-fast-blocks`

## What Changed

JWAS now has an explicit approximate fast-block mode:

```julia
runMCMC(model, phenotypes;
        fast_blocks=[1, 501, 975],
        independent_blocks=true)
```

The default remains:

```julia
independent_blocks=false
```

That default keeps the exact sequential block sweep.
The new independent-block mode freezes the sweep-level corrected phenotype / residual, updates blocks independently, and reconciles all block deltas at the end of the sweep.

`fast_blocks` also accepts explicit marker block starts.
The vector form uses full-sweep chain semantics: one MCMC iteration means one sweep over all supplied blocks.

## Method Coverage

The implementation covers:

| Family | Exact fast blocks | Independent blocks |
| --- | ---: | ---: |
| single-trait BayesA/B/C | yes | yes |
| single-trait BayesR | yes | yes |
| single-trait annotated BayesC | yes | yes |
| single-trait annotated BayesR | yes | yes |
| plain multi-trait BayesA/B/C sampler I | yes | yes |
| plain multi-trait BayesA/B/C sampler II | yes | yes |
| annotated 2-trait BayesC sampler I | yes | yes |
| annotated 2-trait BayesC sampler II | yes | yes |

## Statistical Interpretation

For single-trait models, independent blocks are exact only when:

```text
X_b' W X_c = 0  for all b != c
```

For multi-trait models, residual covariance still couples traits inside each block update.
The independent-block approximation is still about off-block genotype leakage.

So:

- `independent_blocks=false` is the exact default.
- `independent_blocks=true` is an explicit approximation unless the supplied block partition makes off-block weighted genotype crossproducts vanish.
- Pedigree, IBD, recombination, LD, or external block software can help propose better blocks, but the independence assumption is still judged by genotype block leakage.

## Tests Added

The implementation added tests for:

- API plumbing and validation for `independent_blocks`.
- Validation of explicit `fast_blocks` starts.
- Explicit block-start chain-length semantics.
- Single-trait BayesC independent fast-block production path.
- Single-trait annotated BayesC independent fast-block production path.
- Single-trait BayesR independent fast-block production path.
- Single-trait annotated BayesR independent fast-block production path.
- Plain multi-trait BayesC sampler I independent fast-block production path.
- Plain multi-trait BayesC sampler II independent fast-block production path.
- Annotated 2-trait BayesC sampler I independent fast-block production path.
- Annotated 2-trait BayesC sampler II independent fast-block production path.
- A deterministic weighted crossproduct identity test for the exact independent-block condition.
- Two-thread execution smoke coverage for annotated and multi-trait paths.

## Verification Run

Fresh verification completed on this branch:

| Command | Result |
| --- | --- |
| `julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl"); include("test/unit/test_annotated_bayesc.jl"); include("test/unit/test_annotated_bayesr.jl"); include("test/unit/test_multitrait_mcmc.jl")'` | passed, exit code 0 |
| `JULIA_NUM_THREADS=2 julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl"); include("test/unit/test_annotated_bayesr.jl"); include("test/unit/test_multitrait_mcmc.jl")'` | passed, exit code 0 |
| `julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'` | passed, exit code 0 |
| `julia --project=docs --startup-file=no docs/make.jl` | passed, exit code 0 |
| `julia --project=. --startup-file=no test/runtests.jl` | passed, `836/836` tests, `3m11.2s` |

The first full-suite run exposed one stale BayesR test probe that only wrapped the old six-argument `BayesR_block!` signature.
The production call now passes the seventh `independent_blocks` argument, so the probe did not observe dispatch.
The probe was updated to cover both signatures, and the BayesR file plus full suite passed afterward.

## Production Smoke Benchmark

Benchmark driver:

```text
benchmarks/independent_blocks_fast_blocks_smoke.jl
```

Inputs:

- JWAS `demo_7animals` genotype and phenotype fixture.
- `chain_length=20`, `burnin=5`.
- `fast_blocks=[1, 3, 5]`.
- `BLAS.set_num_threads(1)`.
- Production `runMCMC`, not helper-only code.

The smoke benchmark writes:

```text
benchmarks/reports/2026-04-15-independent-blocks-fast-block-parallelism-smoke-results-threads1.csv
benchmarks/reports/2026-04-15-independent-blocks-fast-block-parallelism-smoke-results-threads2.csv
```

### Timing Summary

These numbers verify execution across production paths.
They should not be treated as final speedup estimates because the fixture has only five markers and first-call compilation dominates several exact-mode timings.

| Case | 1 thread exact (s) | 1 thread independent (s) | 2 threads exact (s) | 2 threads independent (s) |
| --- | ---: | ---: | ---: | ---: |
| BayesC | 11.175 | 0.332 | 11.124 | 0.313 |
| BayesR | 0.857 | 0.259 | 1.035 | 0.260 |
| Annotated BayesC | 0.483 | 0.012 | 0.566 | 0.016 |
| Annotated BayesR | 1.737 | 0.240 | 1.792 | 0.246 |
| MT BayesC sampler I | 5.664 | 1.041 | 5.842 | 1.165 |
| MT BayesC sampler II | 0.756 | 0.822 | 0.831 | 1.012 |
| MT annotated BayesC sampler I | 1.183 | 0.788 | 1.319 | 0.918 |
| MT annotated BayesC sampler II | 0.474 | 0.613 | 0.457 | 0.822 |

Interpretation:

- All eight model families ran through both exact fast-block and independent-block production paths.
- The tiny fixture does not show useful two-thread scaling because the threaded work is too small.
- Some independent MT sampler-II smoke runs are slower on this fixture, which is expected when thread overhead exceeds block work.
- A real server benchmark needs larger marker blocks, warmup, multiple seeds, and marker-key-aligned output comparisons.

## Server Usage

Recommended runtime setup:

```bash
export JULIA_NUM_THREADS=<num_cores>
export OPENBLAS_NUM_THREADS=1
julia --project=.
```

Then use:

```julia
runMCMC(model, phenotypes;
        fast_blocks=my_block_starts,
        independent_blocks=true)
```

where `my_block_starts` is a sorted vector beginning at `1`.

## Comparison Policy

This smoke report does not compute marker-effect or posterior-inclusion correlations.
For any future marker-level comparison, marker rows must be aligned by explicit marker ID or SNP index before computing summaries.
Row order or lexicographic sorting alone is not acceptable.

## Practical Conclusion

The implementation now supports the independent-block approximation across the single-trait, annotated single-trait, plain multi-trait, and annotated multi-trait BayesC/BayesR paths covered above.
The next meaningful benchmark should use a larger livestock genotype fixture with user-provided LD or pedigree-informed blocks and report:

- runtime and thread scaling;
- EBV / prediction accuracy;
- marker-effect correlations after marker-key alignment;
- posterior inclusion or any-active recovery after marker-key alignment;
- sensitivity to block partition quality.
