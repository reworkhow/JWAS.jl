# BayesC-only memory implementation

Base: reworkhow master 038c6c8a497ec6e17e957657a5c13ae698377132 (JWAS 2.4.0).
Review branch: performance/bayesc-only-shared-views.

## Included

- Ordinary and fast-block BayesC; single-trait and multi-trait samplers I/II.
- Sequential and independent blocks.
- BayesC-only matrix-adoption defaults and scoped automatic shared storage.
- Contiguous parent-backed training, prediction and block views, with parent-aware
  memory accounting. Other methods retain copying/storage defaults.
- One in-place BayesC block residual path with genotype-typed BLAS multipliers;
  no user residual-mode option and no product workspaces.
- Sequential MT copies only the current alpha block. Residual buffers keep
  their existing element type. Shared ABC samplers select legacy arithmetic
  internally for other methods to keep the scope confined to BayesC.

BayesR and MTBayesR source files, prediction output code, priors, dependencies,
package version and sampling schedules are unchanged from upstream.

## Verification

All 1767 full-suite checks passed, including the unchanged 988 legacy checks, default-versus-copied
storage comparisons, ST/MT residual reconstruction across both precisions,
sequential/independent blocks and three seeds, mixed-precision genotype-typed
multipliers, and preservation of BayesB wrapper arithmetic. Duplicate arithmetic
mode cases were removed with the public switch. No original regression test
assertions were changed or disabled.

All 72 keyed upstream-versus-candidate prediction cases passed absolute/relative
1e-4 tolerance: single-trait plus two-trait samplers I/II, both precisions,
ordinary/sequential/independent updates, three seeds, annotations, partial trait
missingness and held-out animals. Maximum EBV difference was 1.108886e-6;
PEV 8.178929e-7, marker mean 5.1e-7 and marker SD 3.3e-7. Model frequencies
matched exactly. All 18 ordinary cases and 20 of 54 block cases matched exactly.
The unified block arithmetic is not claimed bitwise identical to upstream.

Eight production storage comparisons (one/two traits, ordinary/sequential,
two seeds) reduced retained genotype payload from 230272 to 131584 bytes,
42.9%, when training/prediction outputs overlap. This is not total RSS.

Documentation built successfully from the candidate source; deployment skipped.
Tests used Julia 1.11.7 on Apple Silicon. Optional integration tests and the
full CI runtime/OS matrix were not run. Synthetic prediction checks do not
establish production convergence or equivalence for every dataset.

See benchmarks/reports/2026-10-31-bayesc-only/ for compact keyed evidence and
reproduction details. Existing runtime priors and formulas were preserved.
No master branch or upstream pull request is changed during this review.
