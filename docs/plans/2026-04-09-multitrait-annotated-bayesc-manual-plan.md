# Multi-Trait Annotated BayesC Manual Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a dedicated manual page for dense 2-trait annotated BayesC and link it cleanly from the existing annotated BayesC documentation.

**Architecture:** Keep the current single-trait annotated BayesC page as the general entry point, move multi-trait-specific explanation into a dedicated manual page, and update the manual navigation plus cross-links so users can find the 2-trait workflow directly.

**Tech Stack:** Markdown docs, Documenter.jl navigation, JWAS manual pages

---

### Task 1: Add a dedicated multi-trait manual page

**Files:**
- Create: `docs/src/manual/multitrait_annotated_bayesc.md`

**Step 1: Write the content**

Create a new page that covers:

- dense 2-trait annotated BayesC support status
- four joint states `00`, `10`, `01`, `11`
- 3-step tree prior reconstruction
- startup `Pi` defaults and validation requirements
- `multi_trait_sampler=:auto|:I|:II`
- supported and unsupported runtime modes
- one worked 2-trait example
- output interpretation for joint `pi_<geno>` and annotation step names

**Step 2: Keep examples aligned with production behavior**

Use the documented example style already used in:

- `docs/src/manual/annotated_bayesc.md`
- `docs/src/manual/workflow.md`

Make the example explicit about:

- dense only
- 2 traits only
- joint `Pi` dictionary order by state labels, not by column position

**Step 3: Review for scope**

Confirm the new page does not repeat single-trait block or streaming guidance.

### Task 2: Simplify the existing annotated BayesC page

**Files:**
- Modify: `docs/src/manual/annotated_bayesc.md`

**Step 1: Reduce multi-trait material**

Replace the current long 2-trait sections with:

- a brief note that dense 2-trait annotated BayesC is supported
- a short summary of what is different
- a link to `multitrait_annotated_bayesc.md`

**Step 2: Keep single-trait guidance intact**

Retain:

- single-trait method overview
- single-trait dense example
- `fast_blocks` single-trait guidance
- streaming single-trait guidance

### Task 3: Update manual navigation and cross-links

**Files:**
- Modify: `docs/make.jl`
- Modify: `docs/src/manual/workflow.md`

**Step 1: Add the new page to navigation**

Insert:

- `"Multi-Trait Annotated BayesC" => "manual/multitrait_annotated_bayesc.md"`

under the Manual section in `docs/make.jl`, near the other annotated method guides.

**Step 2: Add workflow link**

Update the method-specific guide list in `docs/src/manual/workflow.md` to include
the new page.

**Step 3: Update workflow text**

Where the workflow currently mentions dense 2-trait annotated BayesC, point
readers to the new dedicated page instead of overloading the single-trait page.

### Task 4: Verify the docs build

**Files:**
- Verify: `docs/src/manual/multitrait_annotated_bayesc.md`
- Verify: `docs/src/manual/annotated_bayesc.md`
- Verify: `docs/src/manual/workflow.md`
- Verify: `docs/make.jl`

**Step 1: Build docs**

Run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected:

- build succeeds
- no broken internal links caused by the new page split

**Step 2: Sanity check the result**

Confirm:

- the new page appears in Manual navigation
- the existing annotated BayesC page links to it
- the workflow page links to it

### Task 5: Record completion

**Files:**
- Create: `docs/plans/2026-04-09-multitrait-annotated-bayesc-manual-implementation.md`

**Step 1: Summarize implementation**

After the docs are updated, write the implementation note capturing:

- files changed
- what moved or was rewritten
- verification command used

**Step 2: Commit**

```bash
git add docs/src/manual/annotated_bayesc.md \
        docs/src/manual/multitrait_annotated_bayesc.md \
        docs/src/manual/workflow.md \
        docs/make.jl \
        docs/plans/2026-04-09-multitrait-annotated-bayesc-manual-design.md \
        docs/plans/2026-04-09-multitrait-annotated-bayesc-manual-plan.md \
        docs/plans/2026-04-09-multitrait-annotated-bayesc-manual-implementation.md
git commit -m "docs: split multi-trait annotated BayesC manual"
```
