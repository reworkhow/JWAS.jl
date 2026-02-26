# SEM Documentation Page Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a single layered SEM manual page (Beginner -> Reference -> Deep Technical), wire it into docs navigation, and cross-link workflow guidance so SEM docs are centralized and up to date.

**Architecture:** Introduce one new docs page (`manual/sem.md`) and structure it in progressive sections with an at-top quick index. Keep workflow docs concise by linking to this canonical SEM page instead of duplicating long SEM guidance. Validate through targeted content checks and full docs build.

**Tech Stack:** Julia (`Documenter.jl` docs build), Markdown docs under `docs/src`, shell checks (`rg`), git.

---

### Task 1: Add Manual navigation entry for SEM page

**Files:**
- Modify: `docs/make.jl`
- Create later in Task 2: `docs/src/manual/sem.md`

**Step 1: Add SEM entry in manual nav**

Modify manual pages list in `docs/make.jl` by inserting:

```julia
"SEM (Beginner to Advanced)" => "manual/sem.md",
```

Place it after `"Workflow"`.

**Step 2: Run docs build to verify failure before page exists (RED)**

Run:

```bash
julia --project=docs docs/make.jl
```

Expected: FAIL with missing file for `manual/sem.md`.

**Step 3: Commit nav-only change**

```bash
git add docs/make.jl
git commit -m "docs: add SEM manual nav entry"
```

### Task 2: Create SEM page skeleton and make docs build pass

**Files:**
- Create: `docs/src/manual/sem.md`

**Step 1: Create minimal page skeleton**

Add:

```markdown
# SEM: Beginner to Advanced

```@contents
Pages = ["sem.md"]
Depth = 2
```

## Quick Index

- [Beginner Quick Start](#Beginner-Quick-Start)
- [Reference](#Reference)
- [Deep Technical Notes](#Deep-Technical-Notes)
- [Troubleshooting](#Troubleshooting)
```

**Step 2: Run docs build to verify pass (GREEN)**

Run:

```bash
julia --project=docs docs/make.jl
```

Expected: PASS and generated `docs/build/manual/sem/index.html`.

**Step 3: Commit skeleton page**

```bash
git add docs/src/manual/sem.md
git commit -m "docs: add SEM page skeleton"
```

### Task 3: Implement Beginner Quick Start section

**Files:**
- Modify: `docs/src/manual/sem.md`

**Step 1: Verify section absent first (RED check)**

Run:

```bash
rg -n "^## Beginner Quick Start$" docs/src/manual/sem.md
```

Expected: non-zero exit (section not yet present).

**Step 2: Add beginner content**

Add a full `## Beginner Quick Start` section containing:
- plain-language SEM summary
- prerequisites/constraints bullets
- copy-paste runnable example using `causal_structure`
- expected output files summary
- quick sanity checks

Include this exact example matrix:

```julia
my_structure = [0.0 0.0 0.0
                1.0 0.0 0.0
                1.0 0.0 0.0]
```

**Step 3: Verify section exists (GREEN check)**

Run:

```bash
rg -n "^## Beginner Quick Start$" docs/src/manual/sem.md
```

Expected: one match.

**Step 4: Build docs**

```bash
julia --project=docs docs/make.jl
```

Expected: PASS.

**Step 5: Commit beginner section**

```bash
git add docs/src/manual/sem.md
git commit -m "docs: add SEM beginner quick start"
```

### Task 4: Implement Reference section

**Files:**
- Modify: `docs/src/manual/sem.md`

**Step 1: Verify reference header absent first (RED check)**

Run:

```bash
rg -n "^## Reference$" docs/src/manual/sem.md
```

Expected: non-zero exit if not added yet.

**Step 2: Add reference tables**

Add `## Reference` with tables for:
- `runMCMC(...; causal_structure=...)` behavior
- validation rules and runtime side effects
- output files (`structure_coefficient_MCMC_samples.txt`, indirect/overall outputs)
- compatibility/limitations

**Step 3: Verify section and required keywords (GREEN checks)**

Run:

```bash
rg -n "^## Reference$|causal_structure|structure_coefficient_MCMC_samples" docs/src/manual/sem.md
```

Expected: matches for header and key terms.

**Step 4: Build docs**

```bash
julia --project=docs docs/make.jl
```

Expected: PASS.

**Step 5: Commit reference section**

```bash
git add docs/src/manual/sem.md
git commit -m "docs: add SEM reference tables"
```

### Task 5: Implement Deep Technical Notes and Troubleshooting

**Files:**
- Modify: `docs/src/manual/sem.md`

**Step 1: Verify deep section absent first (RED check)**

Run:

```bash
rg -n "^## Deep Technical Notes$" docs/src/manual/sem.md
```

Expected: non-zero exit if not added yet.

**Step 2: Add deep technical section**

Add `## Deep Technical Notes` covering:
- notation and matrix construction (`Y`, `Λ`, `Λy`)
- MCMC sampling flow and update points
- indirect/overall marker effect generation pipeline
- caveats/assumptions and source-file pointers

**Step 3: Add troubleshooting section**

Add `## Troubleshooting` with common SEM failures and fixes:
- non-lower-triangular structure
- single-trait misuse
- output filename mismatches tied to genotype-term names
- note that issue #162 regression is covered by unit test

**Step 4: Verify headers and key terms (GREEN checks)**

Run:

```bash
rg -n "^## Deep Technical Notes$|^## Troubleshooting$|issue #162|lower triangular" docs/src/manual/sem.md
```

Expected: matches for both headers and keywords.

**Step 5: Build docs**

```bash
julia --project=docs docs/make.jl
```

Expected: PASS.

**Step 6: Commit deep + troubleshooting sections**

```bash
git add docs/src/manual/sem.md
git commit -m "docs: add SEM deep technical notes and troubleshooting"
```

### Task 6: Add workflow cross-link and finalize

**Files:**
- Modify: `docs/src/manual/workflow.md`
- Verify built output: `docs/build/manual/workflow/index.html` and `docs/build/manual/sem/index.html`

**Step 1: Add concise SEM pointer in workflow page**

In existing SEM mention area, add a short pointer to the new SEM page:

```markdown
For complete SEM usage, reference, and implementation notes, see [SEM: Beginner to Advanced](sem.md).
```

Keep workflow concise and avoid duplicating large SEM content blocks.

**Step 2: Verify cross-link text exists**

Run:

```bash
rg -n "SEM: Beginner to Advanced|sem\.md" docs/src/manual/workflow.md
```

Expected: one or more matches.

**Step 3: Final docs build and smoke checks**

Run:

```bash
julia --project=docs docs/make.jl
rg -n "SEM: Beginner to Advanced" docs/build/manual/workflow/index.html docs/build/manual/sem/index.html
```

Expected: build PASS and both files contain SEM page text/link.

**Step 4: Final commit**

```bash
git add docs/src/manual/sem.md docs/src/manual/workflow.md docs/make.jl
git commit -m "docs: add comprehensive SEM manual page"
```

### Task 7: Optional doc quality pass (if time)

**Files:**
- Modify if needed: `docs/src/manual/sem.md`

**Step 1: Run quick wording/lint pass**

Check for typos and stale wording:

```bash
rg -n "structue|affacts|Catogorical|estimted|computaionally" docs/src/manual/sem.md docs/src/manual/workflow.md
```

Fix only obvious typo-level issues related to SEM docs touched in this plan.

**Step 2: Rebuild docs**

```bash
julia --project=docs docs/make.jl
```

Expected: PASS.

**Step 3: Commit typo-only adjustments**

```bash
git add docs/src/manual/sem.md docs/src/manual/workflow.md
git commit -m "docs: polish SEM manual wording"
```
