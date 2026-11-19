# BayesC no-copy revision implementation

Base: PR #175 head `b8ef24ce2f79e49512891e00f431b4fbc8844a5c`.
Core revision: `2de114ef`.
Design: [BayesC no-copy revision](2026-11-18-bayesc-no-copy-revision-design.md).

## Changes

- Preserved BayesC matrix-adoption and matching-precision shared-view defaults,
  block views, priors, sampling schedules and sampler arithmetic.
- Included initial Matrix precision in automatic sharing eligibility. Default
  Float32 import followed by explicitly requested double-precision MCMC reaches
  the existing promotion/copying path. Explicit incompatible sharing still fails
  before conversion.
- Added an internal conversion-buffer estimate that deduplicates retained dense
  genotype parents across terms and counts each required destination. The existing
  memory guard checks this estimate before Float64 genotype allocation. It is not
  a bound on total process RSS or memory outside the known genotype terms.
- Made subsequent dense marker estimates use the actual stored element type.
- Skipped the genotype matrix slice when every SNP passes QC. Actual marker
  removal and marker/allele-frequency/annotation selection semantics are retained.
- Updated the public manual and added isolated-process production-path RSS
  measurement scripts with strict keyed comparisons and source fingerprints.

## Verification

The new regression file was run before production edits: the genuine RED run
reported 30 passes, 12 failures and 2 errors, exposing the intended regressions
and the missing conversion estimator. After the fixes and final annotation
metadata assertions, the focused run passed all 57 checks.

The complete default test suite passed 1,824/1,824 checks, including the unchanged
legacy coverage. Optional integration tests remained disabled by the repository's
default configuration. The documentation build passed; deployment was skipped.

Commands used from the isolated worktree:

```sh
GKSwstype=100 JULIA_PKG_OFFLINE=true julia --project=. --startup-file=no --compiled-modules=existing test/runtests.jl
GKSwstype=100 JULIA_PKG_OFFLINE=true JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=docs --startup-file=no --compiled-modules=existing docs/make.jl
python3 -B benchmarks/run_bayesc_no_copy_rss.py --self-test
python3 -B benchmarks/run_bayesc_no_copy_rss.py
```

The ignored documentation environment was pointed at this worktree's source;
a writable temporary first Julia depot was used for documentation runtime logs.
Julia 1.11.7 on Apple Silicon was used. Dependencies and package version were not
changed. Independent scope and code-quality reviews approved the core changes.

## Benchmark interpretation and remaining work

The expanded report is under
`benchmarks/reports/2026-11-18-bayesc-no-copy-revision/`. It records separate
genotype-parent payloads, caller-held buffers and process-lifetime peak RSS,
including package loading, compilation and fixture preparation. Its source
fingerprint identifies the tested production code and benchmark runners.

All 36 fresh-process cases completed, and all 24 strictly keyed comparisons
matched exactly for exported EBV, PEV, marker mean, SD and inclusion frequency.
Joined animal and marker samples were inspected after key/count validation.
Every adopted case retained the caller buffer, including QC-on cases, and kept
one model genotype parent. At the largest measured size (2,048 by 4,097), model
genotype parents occupied 33,562,624 bytes in the adopted/shared path versus
58,734,592 bytes in the copied path. The report separately accounts for retained
caller storage and measured process peaks; these payload figures are not RSS.

Practical increasing synthetic sizes and repeated seeds are used; 50,000
individuals by millions of SNPs were not executed. Arithmetic storage-size
examples are extrapolations, not measured RSS or demonstrated production
feasibility. Cold-process RSS measurements do not establish a speedup or a
reduction in total peak memory for every workload.

CSV/DataFrame loading still creates intermediate representations, and actual
QC filtering still allocates a replacement genotype matrix. The earlier import
can therefore exhaust memory before either runMCMC guard executes. These larger
loading/filtering changes, new packed/streaming functionality, and any redesign
of storage versus computation precision remain separate work.
