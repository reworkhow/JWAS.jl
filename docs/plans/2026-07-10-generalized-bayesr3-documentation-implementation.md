# Generalized BayesR3 Documentation Implementation

Date: 2026-07-11

## Scope and Terminology Decisions

This documentation-only change reframes `fast_blocks` as JWAS's generalized
BayesR3 block strategy. The strategy consists of block right-hand-side
construction, repeated conditional marker updates through the block Gram
matrix, and block-exit residual reconciliation. The documentation credits the
original BayesR3 strategy without claiming that JWAS reproduces every model or
scheduling detail of the original software.

The revised manual treats the following as independent choices:

1. statistical model;
2. block partition;
3. repeat policy; and
4. execution schedule.

It distinguishes the exact sequential block schedule from the explicitly
approximate `independent_blocks=true` schedule, except when off-block weighted
genotype crossproducts vanish. It also separates JWAS's arithmetic
operation-count model from the BayesR3 paper's empirical runtime fit; the two
expressions are not presented as a software-to-software benchmark.

The repeat-count comparison now uses `b` for nominal block size, `s_i` for the
realized size of block `i`, and `r_i` for its inner repeat count. It records the
different repeat policies used by the BayesR3 paper and the current JWAS BayesC
and BayesR implementations without treating those scheduling differences as a
different block linear-algebra strategy.

## Files Changed

- `docs/src/manual/block_bayesc.md`: retitled and reorganized the page around
  the BayesR3 block strategy; documented supported model paths, independent
  design dimensions, explicit unequal blocks, repeat policies, exact versus
  independent-block schedules, and the distinction between operation counts
  and the published empirical runtime fit. The existing path was retained so
  inbound links remain valid.
- `docs/src/manual/annotated_bayesr.md`: aligned the annotated BayesR discussion
  with the generalized strategy and exact transition-schedule terminology.
- `docs/src/manual/annotated_bayesc.md`,
  `docs/src/manual/multitrait_annotated_bayesc.md`,
  `docs/src/manual/bayesc_bayesr_comparison.md`, and
  `docs/src/manual/sem.md`: updated visible cross-reference labels to
  `BayesR3 Block Strategy` while preserving their link target.
- `docs/src/manual/workflow.md`: broadened `fast_blocks` wording from a
  BayesC-only path to the supported generalized strategy and aligned the
  cross-reference label.
- `docs/make.jl`: renamed the navigation entry to
  `BayesR3 Block Strategy` while retaining `manual/block_bayesc.md` as its
  destination.
- `docs/plans/2026-07-10-generalized-bayesr3-documentation-implementation.md`:
  recorded the completed documentation scope, terminology decisions,
  verification results, and remaining explicit-block limitation.

## Behavior and API

Sampler behavior and the public API are unchanged. No production sampler/source
files under `src/` or test files were modified; `docs/make.jl` changed only the
navigation label. The existing forms `fast_blocks=true`,
`fast_blocks=<integer>`, `fast_blocks=<ordered start positions>`, and
`independent_blocks=true` retain their current behavior.

## Verification

The final documentation branch was checked against base commit `48de5e84`.

- `git diff --check` completed successfully with no whitespace errors.
- `rg -n -i "JWAS uses BayesC \\(not BayesR\\)|accelerated approximation|BayesR3 paper fit vs standard|apparent.*gap" docs/src docs/make.jl`
  returned no matches.
- `/Users/haocheng/.juliaup/bin/julia --project=docs --startup-file=no docs/make.jl`
  completed successfully. Documenter built cross-references and rendered the
  HTML site without missing-page or unresolved-link errors. Its only warning
  was that it could not detect a deployment environment and therefore skipped
  deployment, which is expected for this local build. The explicit Julia path
  was used because `julia` was not available on the shell `PATH`.
- Rendered-output searches and direct inspection of
  `docs/build/manual/block_bayesc/index.html` confirmed the page title and
  navigation label, user-defined unequal-block guidance, the empirical runtime
  fit section, and the `independent_blocks` distinction. Inspection of rendered
  cross-page links confirmed that `BayesR3 Block Strategy` resolves to the
  retained `block_bayesc` page.
- `docs/build/` remained ignored and was not added to the branch.

## Remaining Limitation

The current explicit-block API accepts ordered block start positions. Its
blocks are therefore contiguous in the current marker order, although they may
be unequal in size. Arbitrary non-contiguous marker membership is not currently
supported.
