# BayesC no-copy revision design

Base: PR #175 head `b8ef24ce2f79e49512891e00f431b4fbc8844a5c`.

## Objective and scope

Preserve BayesC matrix adoption, contiguous training/prediction views, marker-block
views and the current sampler arithmetic. Fix precision eligibility and check
known conversion-buffer memory before allocation. Avoid QC's genotype copy when
all SNPs pass. Keep CSV ingestion redesign, memory-efficient actual SNP filtering,
new storage formats and a numerical precision API redesign separate.

For 50,000 individuals and one million SNPs, a Float32 genotype matrix is 200 GB
decimal. Float64 promotion requires a 400 GB destination while the source remains
alive. Checks must account for actual buffer types and report their limits; they
are not a guarantee about total process RSS.

## Behavior

- Keep `copy_genotypes=false` for BayesC matrix inputs and the existing automatic
  sharing default for matching-precision supported models.
- Automatic sharing additionally requires the initial Matrix precision to match
  the requested run precision. A mismatched default import follows the existing
  requested conversion/copying path. Explicit `shared_genotypes=true` retains its
  early mismatch validation. No new public conversion option is introduced.
- Before an existing Float64 genotype promotion, reuse `memory_guard` and
  `memory_guard_ratio` to check unique retained dense genotype parents plus the
  new destination buffer. A failing `:error` check must leave the genotype buffer
  unconverted. `:warn` and `:off` retain their current meanings.
- Base dense memory accounting on the actual stored element type.
- When QC selects every SNP, keep the genotype matrix itself. Continue the same
  marker/allele-frequency/annotation selection semantics and retain the existing
  allocating behavior when SNPs are actually removed.

## Verification and evidence

Write failing focused tests first, observe the precision regression, missing
early guard, and all-pass QC alias failure, then implement the fixes. Verify
parent identity rather than brittle exact allocation counts. Retain existing
keyed numerical comparisons and multi-seed coverage. Run the full test suite and
docs build against this worktree's JWAS source.

Extend the production-path benchmark with isolated-process peak RSS measurements
through import, QC, alignment and MCMC at practical increasing sizes. Include QC
on/off and copied/shared matrix paths plus the unchanged file loader. Report
retained buffers separately, use repeated seeds, validate marker keys/counts and
joined samples before numerical comparisons, and label scale extrapolations.

## Implementation sequence

1. Add and run regression tests; implement sharing eligibility, early conversion
   guard and all-pass QC shortcut. Run focused tests.
2. Extend and execute the benchmark; update the public manual and compact report.
3. Independently review code and scope, run all tests and docs, record actual
   results in an implementation note, commit and update the existing PR branch.
