# JWAS Agent Notes

## Feature Workflow
- For new feature work, create `docs/plans/YYYY-MM-DD-<feature>-design.md` before implementation.
- After implementation, create `docs/plans/YYYY-MM-DD-<feature>-implementation.md`.
- Keep both files in git as part of the feature record.

## Code Changes
- Prefer small, reviewable commits.
- Preserve existing JWAS behavior unless the task explicitly changes it.
- For new methods or major behavior changes, add or update tests before merge.

## Verification
- Run the relevant test coverage before claiming completion.
- For larger changes, run `julia --project=. --startup-file=no test/runtests.jl`.
- If docs are changed, also run `julia --project=docs --startup-file=no docs/make.jl`.

## Benchmarking
- Benchmark the production JWAS path, not only prototype or helper code.
- For stochastic methods, do not rely on a single short chain; use longer runs or multiple seeds.
- Save benchmark reports under `benchmarks/reports/` and design/implementation notes under `docs/plans/`.

### Result Comparison
- Treat benchmark comparisons as production work. A wrong comparison is worse than no comparison.
- Never compare marker-level outputs by row order or lexicographic sort alone. Always join by an explicit key such as `Marker_ID` or SNP index.
- Normalize marker IDs before joining. Convert quoted IDs, numeric indices, and string labels into one canonical key first.
- Before computing any correlation or summary, assert all of the following:
  - both sides have the same expected marker count
  - the marker-key sets match
  - the join produced the expected row count
  - no duplicated or dropped markers remain
- Save comparison-ready tables with explicit keys and aligned values when possible.
- Before reporting a result, inspect a small sample of joined rows to confirm the alignment is real.
