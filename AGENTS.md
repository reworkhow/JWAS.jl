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
