# Release Polish Documentation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Polish JWAS documentation so recent annotated, streaming, and fast-block updates are discoverable and described consistently.

**Architecture:** This is a docs-only pass. Update release notes and entry points first, then align manual pages around exact dense/block paths, approximate independent-block paths, annotated method boundaries, and streaming constraints. Do not change Julia source, benchmark scripts, or generated assets.

**Tech Stack:** Julia Documenter docs in `docs/`, Markdown manual pages, package `CHANGELOG.md`, package `README.md`.

---

### Task 1: Update Release Notes

**Files:**
- Modify: `CHANGELOG.md`

**Step 1: Read recent implementation notes**

Run:

```bash
sed -n '1,220p' docs/plans/2026-04-15-independent-blocks-fast-block-parallelism-implementation.md
sed -n '1,180p' docs/plans/2026-04-09-multitrait-annotated-bayesc-manual-implementation.md
sed -n '1,180p' docs/plans/2026-04-02-simulated-annotations-dataset-implementation.md
```

Expected: identify recent user-facing changes without adding benchmark claims not present in the docs.

**Step 2: Add a current changelog section**

Add a new top section to `CHANGELOG.md` above `2026-02-26`.

Use this structure:

```markdown
## 2026-05-01

### Added

- Documented independent-block fast-block sampling for supported dense BayesC, BayesR, annotated BayesC, annotated BayesR, and multi-trait BayesC paths.
- Documented dense 2-trait annotated BayesC as a first-class workflow with explicit sampler-selection guidance.
- Improved documentation for simulated annotation datasets and streaming genotype workflows.

### Documentation

- Clarified exact sequential fast blocks versus approximate `independent_blocks=true` updates.
- Clarified support boundaries for annotated methods and streaming genotype storage.
```

Adjust wording based on the implementation notes, but keep it concise.

**Step 3: Commit checkpoint**

Do not commit yet if subsequent docs are being edited in the same small batch. Commit after Task 2 unless this task becomes large.

---

### Task 2: Refresh README And Docs Landing Page

**Files:**
- Modify: `README.md`
- Modify: `docs/src/index.md`

**Step 1: Update README positioning**

In `README.md`, replace the short one-sentence package description with a concise paragraph that mentions:

- Bayesian mixed models for genomic prediction and GWAS
- single-trait and multi-trait analyses
- dense, streaming, and fast-block genotype workflows
- annotated BayesC/BayesR support

Keep installation and author links intact.

**Step 2: De-emphasize old wiki examples**

In `README.md`, change the examples bullet from wiki-first wording to docs-first wording. Keep the wiki link as legacy/community examples, but direct current users to the Documenter manual.

**Step 3: Update docs landing feature list**

In `docs/src/index.md`, add bullets for:

- annotated BayesC and BayesR
- dense 2-trait annotated BayesC
- streaming genotype storage for memory-sensitive runs
- exact fast blocks and approximate independent blocks

Fix obvious stale wording while staying concise. If you touch the "singe-step" typo, do not add single-step debugging content.

**Step 4: Run text sanity check**

Run:

```bash
rg -n "singe-step|wiki|independent_blocks|Annotated|streaming|fast_blocks" README.md docs/src/index.md
```

Expected: no accidental stale "singe-step"; wiki wording should no longer be the primary route for examples.

**Step 5: Commit checkpoint**

Run:

```bash
git add CHANGELOG.md README.md docs/src/index.md
git commit -m "docs: refresh release overview"
```

---

### Task 3: Align Annotated Method Manual Pages

**Files:**
- Modify: `docs/src/manual/annotated_bayesc.md`
- Modify: `docs/src/manual/annotated_bayesr.md`
- Modify: `docs/src/manual/multitrait_annotated_bayesc.md`
- Modify: `docs/src/manual/public.md`

**Step 1: Clarify single-trait annotated BayesC page**

In `docs/src/manual/annotated_bayesc.md`:

- keep the page single-trait first
- ensure it links to `multitrait_annotated_bayesc.md` for dense 2-trait workflows
- mention dense, fast-block, independent-block, and streaming support only where accurate
- avoid duplicating the multi-trait annotated BayesC manual

**Step 2: Clarify annotated BayesR page**

In `docs/src/manual/annotated_bayesr.md`:

- ensure supported storage/block modes are stated consistently with the current implementation
- align wording with annotated BayesC for annotation matrix validation and prior interpretation
- avoid implying unsupported multi-trait annotated BayesR support

**Step 3: Clarify multi-trait annotated BayesC page**

In `docs/src/manual/multitrait_annotated_bayesc.md`:

- keep support boundaries visible near the top
- state dense 2-trait support explicitly
- state `multi_trait_sampler=:auto`, `:I`, and `:II` choices
- mention fast-block and independent-block support only if the page already documents the current implemented coverage

**Step 4: Check public interface page**

In `docs/src/manual/public.md`, verify current public functions include the key user entry points:

- `get_genotypes`
- `prepare_streaming_genotypes`
- `runMCMC`
- `outputEBV`
- `outputMCMCsamples`

Add missing public docs entries only if the docstring exists and Documenter can build it.

**Step 5: Search for contradictions**

Run:

```bash
rg -n "multi_trait_sampler|independent_blocks|fast_blocks|storage=:stream|Annotated BayesR|Annotated BayesC" docs/src/manual/annotated_bayesc.md docs/src/manual/annotated_bayesr.md docs/src/manual/multitrait_annotated_bayesc.md docs/src/manual/public.md
```

Expected: repeated claims are consistent across pages.

**Step 6: Commit checkpoint**

Run:

```bash
git add docs/src/manual/annotated_bayesc.md docs/src/manual/annotated_bayesr.md docs/src/manual/multitrait_annotated_bayesc.md docs/src/manual/public.md
git commit -m "docs: align annotated method guidance"
```

---

### Task 4: Clarify Fast-Block And Streaming Guidance

**Files:**
- Modify: `docs/src/manual/block_bayesc.md`
- Modify if needed: `docs/src/manual/streaming_genotype_backend.md`
- Modify if needed: `docs/src/manual/large_genotype_data_streaming.md`
- Modify if needed: `docs/src/manual/streaming_genotype_walkthrough.md`

**Step 1: Tighten block BayesC guidance**

In `docs/src/manual/block_bayesc.md`, ensure these points are easy to find:

- `fast_blocks=true` or numeric fast blocks keep exact sequential block behavior when `independent_blocks=false`
- explicit block starts are marker start positions and must begin at `1`
- explicit block starts use full-sweep chain semantics
- `independent_blocks=true` is approximate unless off-block weighted genotype crossproducts are effectively zero
- server users should set `OPENBLAS_NUM_THREADS=1` when using Julia threads

**Step 2: Check streaming pages for consistency**

Search streaming docs:

```bash
rg -n "storage=:stream|double_precision|single_step_analysis|multi-trait|fast_blocks|BayesC" docs/src/manual/streaming_genotype_backend.md docs/src/manual/large_genotype_data_streaming.md docs/src/manual/streaming_genotype_walkthrough.md
```

Only edit these pages if they contain stale or contradictory wording. Keep edits small.

**Step 3: Commit checkpoint**

Run:

```bash
git add docs/src/manual/block_bayesc.md docs/src/manual/streaming_genotype_backend.md docs/src/manual/large_genotype_data_streaming.md docs/src/manual/streaming_genotype_walkthrough.md
git commit -m "docs: clarify block and streaming workflows"
```

If no streaming pages changed, omit them from `git add`.

---

### Task 5: Build Docs And Write Implementation Note

**Files:**
- Create: `docs/plans/2026-05-01-release-polish-implementation.md`

**Step 1: Build docs**

Run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected: exit code 0. If Documenter reports missing docstrings from `public.md`, either add the docstring in a separate scoped code-doc commit or remove the new `@docs` entry.

**Step 2: Final search**

Run:

```bash
rg -n "singe-step|old wiki|TODO|FIXME|independent block|independent_blocks|storage=:stream" README.md CHANGELOG.md docs/src
```

Expected: no accidental stale wording in changed files. Existing unrelated TODO/FIXME hits can be left alone if they are not from this task.

**Step 3: Write implementation note**

Create `docs/plans/2026-05-01-release-polish-implementation.md` with:

- summary of changed docs
- key wording decisions
- verification command and result
- note that production benchmarks were not rerun
- note that unrelated untracked benchmark/debug artifacts were not touched

**Step 4: Commit final docs**

Run:

```bash
git add docs/plans/2026-05-01-release-polish-implementation.md
git commit -m "docs: record release polish implementation"
```

If Task 5 includes fixes to docs from the build step, include those exact files in the commit.

---

### Task 6: Final Branch Review

**Files:**
- Inspect all changed files.

**Step 1: Show branch commits**

Run:

```bash
git log --oneline --decorate origin/master..HEAD
```

Expected: small docs-focused commits.

**Step 2: Show final status**

Run:

```bash
git status -sb
```

Expected: no tracked modifications. Existing unrelated untracked benchmark/debug artifacts may still be present and should remain uncommitted.

**Step 3: Summarize**

Report:

- files changed
- docs build result
- commits created
- any skipped verification
