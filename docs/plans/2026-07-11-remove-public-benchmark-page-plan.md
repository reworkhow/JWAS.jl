# Remove Public Benchmark Page Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove the extrapolated fast-block benchmark from the public JWAS manual while preserving its measured calibration record under `benchmarks/reports/` with accurate limitations.

**Architecture:** Move the existing Markdown content out of `docs/src/manual/` so the measurement record remains versioned but is no longer presented as public performance guidance. Remove all Documenter navigation, index, and manual links to the page, then verify that the public site builds without the benchmark page and that the internal report clearly distinguishes measured calibration timings from extrapolated target estimates.

**Tech Stack:** Markdown, Documenter.jl, Git, Julia 1.11 documentation environment.

---

### Task 1: Preserve the benchmark as an internal calibration report

**Files:**
- Move: `docs/src/manual/benchmark.md` to `benchmarks/reports/2026-02-26-fast-block-calibration-and-target-extrapolation.md`

**Step 1: Inspect the source page**

Run:

```bash
sed -n '1,180p' docs/src/manual/benchmark.md
```

Expected: the page calls itself a real target-scale benchmark while later notes reveal that target times were extrapolated from smaller calibration runs.

**Step 2: Move the file with Git history preserved**

Run:

```bash
git mv docs/src/manual/benchmark.md \
  benchmarks/reports/2026-02-26-fast-block-calibration-and-target-extrapolation.md
```

**Step 3: Rewrite the report framing**

Use the heading:

```markdown
# Fast-Block Calibration Timings and Target Extrapolation
```

The opening must state that:

- wall-clock timings at `P=100,000` and `P=200,000` were measured;
- the `P=2,000,000`, `chain_length=2000` values and speed ratio were extrapolated;
- this is not a direct full-scale benchmark or a general performance guarantee;
- the two Slurm jobs ran on different CPU node types.

Keep the measured tables, extrapolation procedure, and results intact as an engineering record. Rename “Summary Comparison” to “Extrapolated Target Comparison” and qualify the `213x` value as an extrapolated ratio for this calibration setup.

**Step 4: Verify the internal report language**

Run:

```bash
rg -n -i "measured|extrapolat|direct full-scale|guarantee|different CPU" \
  benchmarks/reports/2026-02-26-fast-block-calibration-and-target-extrapolation.md
git diff --check
```

Expected: measured and extrapolated evidence are unambiguous and there are no whitespace errors.

**Step 5: Commit**

```bash
git add docs/src/manual/benchmark.md \
  benchmarks/reports/2026-02-26-fast-block-calibration-and-target-extrapolation.md
git commit -m "docs: archive fast-block calibration benchmark"
```

### Task 2: Remove the benchmark from public documentation

**Files:**
- Modify: `docs/make.jl:25`
- Modify: `docs/src/index.md:70`
- Modify: `docs/src/manual/block_bayesc.md:9`
- Modify: `docs/src/manual/workflow.md:318`

**Step 1: Remove public navigation and index entries**

Delete:

```julia
"Benchmark" => "manual/benchmark.md",
```

from `docs/make.jl`, and remove `"manual/benchmark.md"` from the source-page list in `docs/src/index.md`.

**Step 2: Remove public links**

Remove the sentences linking to the benchmark from the BayesR3 block-strategy page and workflow page. Do not replace them with another speed claim.

Retain existing practical advice that users should benchmark their own block sizes and workloads.

**Step 3: Search for stale public references**

Run:

```bash
if rg -n "benchmark\.md|\[Benchmark\]|real cluster benchmark|target-scale benchmark" \
  docs/make.jl docs/src; then
  exit 1
fi
```

Expected: no public documentation or navigation references remain.

**Step 4: Inspect scope and formatting**

Run:

```bash
git diff --check
git diff -- docs/make.jl docs/src
```

Expected: only the intended navigation, index, and two link sentences changed.

**Step 5: Commit**

```bash
git add docs/make.jl docs/src/index.md \
  docs/src/manual/block_bayesc.md docs/src/manual/workflow.md
git commit -m "docs: remove extrapolated benchmark from manual"
```

### Task 3: Verify and record the documentation correction

**Files:**
- Create: `docs/plans/2026-07-11-remove-public-benchmark-page-implementation.md`
- Verify: `docs/build/`

**Step 1: Run final source checks**

Run:

```bash
git diff --check HEAD~2..HEAD
test -f benchmarks/reports/2026-02-26-fast-block-calibration-and-target-extrapolation.md
test ! -f docs/src/manual/benchmark.md
if rg -n "benchmark\.md|\[Benchmark\]|real cluster benchmark|target-scale benchmark" \
  docs/make.jl docs/src; then
  exit 1
fi
```

Expected: the internal report exists, the public page is absent, and no stale public links remain.

**Step 2: Build the documentation**

Run:

```bash
/Users/haocheng/.juliaup/bin/julia --project=docs --startup-file=no docs/make.jl
```

Expected: Documenter exits successfully; the local deployment-skipped warning is acceptable.

**Step 3: Inspect rendered navigation**

Run:

```bash
if rg -n "manual/benchmark|>Benchmark<|target-scale benchmark" docs/build; then
  exit 1
fi
```

Expected: the generated public site contains no benchmark page or navigation label.

**Step 4: Write the implementation record**

Create `docs/plans/2026-07-11-remove-public-benchmark-page-implementation.md` documenting:

- why the page was removed from public guidance;
- the internal report destination;
- files changed;
- confirmation that no Julia source, API, sampler, test, or benchmark data changed;
- verification commands and outcomes;
- the distinction between measured calibration timings and extrapolated target estimates.

**Step 5: Review final scope**

Run:

```bash
git status --short
git diff --stat 6b1cf84e..HEAD
git diff 6b1cf84e..HEAD -- docs/make.jl docs/src benchmarks/reports \
  docs/plans/2026-07-11-remove-public-benchmark-page-implementation.md
```

Expected: only the approved documentation correction and records are present; ignored build output and unrelated user files are absent.

**Step 6: Commit**

```bash
git add docs/plans/2026-07-11-remove-public-benchmark-page-implementation.md
git commit -m "docs: record public benchmark page removal"
```
