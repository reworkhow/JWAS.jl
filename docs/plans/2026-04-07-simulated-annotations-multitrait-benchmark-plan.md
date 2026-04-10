# Simulated Annotations Multi-Trait Benchmark Plan

**Goal:** Add a native 2-trait simulated annotation fixture and benchmark dense
annotated 2-trait BayesC against the requested baselines on the packaged
dataset.

**Architecture:** Keep the existing single-trait fixture files untouched, add
parallel 2-trait files in the same packaged dataset directory, then build a
benchmark harness that resolves those files through `JWAS.Datasets`, runs the
production methods, joins outputs by explicit marker IDs, and writes both raw
and aggregated summaries.

**Tech Stack:** R for fixture regeneration, Julia for JWAS benchmarks and
reporting, CSV/DataFrames/JWAS.Datasets.

---

### Task 1: Extend the packaged fixture with native 2-trait files

**Files:**
- Modify: `src/4.Datasets/data/simulated_annotations/generate_dataset.R`
- Modify: `src/4.Datasets/data/simulated_annotations/README.md`
- Create/update generated files in `src/4.Datasets/data/simulated_annotations/`

**Steps:**

1. Keep the current single-trait files unchanged.
2. Add 2-trait simulation logic that writes:
   - `phenotypes_mt.csv`
   - `annotations_mt.csv`
   - `truth_mt.csv`
3. Document the state allocation, effect generation, and annotation generation
   in both the generator and README.

### Task 2: Add regression coverage for the new fixture files

**Files:**
- Modify: `test/unit/test_misc_coverage.jl`

**Steps:**

1. Verify `Datasets.dataset(...)` resolves the new 2-trait files.
2. Optionally sanity-check expected columns to protect the benchmark input
   contract.

### Task 3: Implement the multi-trait simulated benchmark harness

**Files:**
- Create: `benchmarks/simulated_annotations_multitrait_comparison.jl`

**Steps:**

1. Resolve the packaged 2-trait files through `JWAS.Datasets`.
2. Implement the requested method matrix:
   - dense multi-trait BayesC
   - dense multi-trait annotated BayesC
   - dense single-trait BayesC / annotated BayesC on each trait
   - dense BayesR on each trait
   - optionally dense annotated BayesR on each trait as an additional baseline
3. For each run, save:
   - runtime
   - trait-wise prediction correlation
   - trait-wise effect correlation
   - trait-wise top-k recall
   - any-active recall where applicable
   - annotation coefficient summaries for annotated methods
4. Write raw and aggregate CSV summaries.

### Task 4: Run the benchmark and write the report

**Files:**
- Create: `benchmarks/reports/2026-04-07-simulated-annotations-multitrait-report.md`
- Create: `docs/plans/2026-04-07-simulated-annotations-multitrait-benchmark-implementation.md`

**Steps:**

1. Run the benchmark with multiple seeds and production JWAS paths only.
2. Save the raw outputs under a benchmark output directory.
3. Write the markdown report from the actual saved summaries.
4. Record the implementation details and commands in the implementation note.

### Task 5: Verification

**Commands:**

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
julia --project=. --startup-file=no benchmarks/simulated_annotations_multitrait_comparison.jl <OUTPUT_DIR>
julia --project=docs --startup-file=no docs/make.jl
```
