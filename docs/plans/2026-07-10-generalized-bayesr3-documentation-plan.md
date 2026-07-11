# Generalized BayesR3 Documentation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reframe the JWAS fast-block manual as a generalized BayesR3 strategy that supports multiple Bayesian models and user-defined unequal contiguous blocks while preserving accurate scheduling, exactness, and performance language.

**Architecture:** Keep `docs/src/manual/block_bayesc.md` at its existing path so inbound links remain valid, but broaden its heading and organization from a BayesC-only page to a strategy-level page. Align navigation and cross-references, make model, partition, repeat count, and execution schedule independent concepts, and verify the rendered Documenter site without changing source behavior.

**Tech Stack:** Markdown, Documenter.jl, Julia 1.10 project environments, repository documentation build.

---

### Task 1: Reframe the main page around the generalized BayesR3 strategy

**Files:**
- Modify: `docs/src/manual/block_bayesc.md:1-98`

**Step 1: Capture the stale terminology before editing**

Run:

```bash
rg -n "Block BayesC|similar to the blocked update|When This Path Is Used|Single-Trait Block BayesC" docs/src/manual/block_bayesc.md
```

Expected: the page is titled `Block BayesC`, describes the implementation as merely similar to BayesR3, and organizes the introduction around BayesC.

**Step 2: Replace the heading and introductory contract**

Use this heading:

```markdown
# BayesR3 Block Strategy in JWAS (`fast_blocks`)
```

State explicitly that:

- JWAS adopts and generalizes the BayesR3 computational strategy;
- the strategy is block RHS construction, repeated conditional marker updates through the block Gram matrix, and block-exit residual reconciliation;
- JWAS does not claim to reproduce every model or scheduling detail of the original BayesR3 software;
- the strategy is used by supported BayesA/B/C, BayesR, annotated, and multi-trait paths listed on the page.

**Step 3: Add the independent design dimensions**

Add a compact section that separately defines:

```markdown
1. Statistical model
2. Block partition
3. Repeat policy
4. Execution schedule
```

Explain that changing the model from BayesR to BayesC does not change the core block linear-algebra strategy, and that changing the block partition does not by itself select the independent-block approximation.

**Step 4: Preserve and clarify the existing API examples**

Keep examples for:

```julia
fast_blocks=true
fast_blocks=64
fast_blocks=[1, 501, 975, 1600]
independent_blocks=true
```

Clarify that an explicit vector contains ordered marker start positions and therefore defines unequal contiguous blocks in current marker order.

**Step 5: Review the first section**

Run:

```bash
sed -n '1,150p' docs/src/manual/block_bayesc.md
```

Expected: a reader can identify the strategy, the supported model families, and the three input forms without inferring that JWAS is BayesC-only or an exact reproduction of BayesR3.

**Step 6: Commit the main-page framing**

```bash
git add docs/src/manual/block_bayesc.md
git commit -m "docs: frame fast blocks as generalized BayesR3 strategy"
```

### Task 2: Document user-defined unequal blocks and convergence implications

**Files:**
- Modify: `docs/src/manual/block_bayesc.md:26-95`

**Step 1: Add a user-defined block section**

Document that users can choose boundaries from LD, recombination, IBD, pedigree information, or genomic regions and that blocks need not have equal sizes.

State the current constraint precisely:

```markdown
The current API accepts ordered block start positions, so every block is contiguous in the current marker order. Arbitrary non-contiguous marker membership is not currently supported.
```

**Step 2: Qualify the mixing statement**

Explain that grouping strongly coupled markers can improve locality or mixing per unit of computation, but do not promise improved convergence. Name the relevant dependencies: posterior correlations, block sizes, boundaries, and repeat policy.

**Step 3: Check for overclaims**

Run:

```bash
rg -n -i "guarantee|improves convergence|faster convergence|user-defined|unequal|contiguous" docs/src/manual/block_bayesc.md
```

Expected: unequal and contiguous semantics are explicit, and no sentence guarantees improved convergence.

**Step 4: Commit the user-block guidance**

```bash
git add docs/src/manual/block_bayesc.md
git commit -m "docs: explain informed unequal marker blocks"
```

### Task 3: Correct the BayesR3 comparison and repeat-count notation

**Files:**
- Modify: `docs/src/manual/block_bayesc.md:172-219`
- Reference: `src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl:145-155`
- Reference: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl:22-25`

**Step 1: Replace the broad BayesC statement**

Replace:

```markdown
JWAS uses BayesC (not BayesR)
```

with wording that identifies the table as a worked comparison of JWAS block BayesC while noting that JWAS also applies the strategy to BayesR and supported annotated paths.

**Step 2: Define block and repeat notation**

Use:

```markdown
- `b`: nominal block size
- `s_i`: realized size of block `i`
- `r_i`: inner repeat count for block `i`
```

**Step 3: State the original BayesR3 convention**

Document that most BayesR3 blocks have size `b`, a final block may have `s_i < b`, and BayesR3 still uses `r_i = b` for that block. Prefer “every SNP effect receives the same number of updates under the paper's chain-length convention” over “separate marker chains have the same length.”

**Step 4: State the JWAS conventions**

Document that:

- BayesC uses `r_i = s_i`;
- BayesR uses `r_i = 1` during burn-in and `r_i = s_i` afterward;
- explicit unequal block starts therefore allow different repeat counts among blocks;
- these schedules differ from original BayesR3 without ceasing to use the generalized block strategy.

**Step 5: Verify wording against source**

Run:

```bash
rg -n "nreps|bayesr_block_nreps" \
  src/1.JWAS/src/markers/BayesianAlphabet/BayesABC.jl \
  src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl \
  docs/src/manual/block_bayesc.md
```

Expected: the documentation matches the current BayesC and BayesR repeat logic.

**Step 6: Commit the comparison correction**

```bash
git add docs/src/manual/block_bayesc.md
git commit -m "docs: distinguish BayesR3 and JWAS repeat policies"
```

### Task 4: Separate exact scheduling from explicit approximation

**Files:**
- Modify: `docs/src/manual/block_bayesc.md:95-170`
- Modify: `docs/src/manual/annotated_bayesr.md:346-354`

**Step 1: Tighten the exact sequential definition**

State that sequential blocks compose conditional Gibbs updates and reconcile the global residual before the next block. Clarify that this is a different transition schedule from dense marker-by-marker Gibbs, but is the exact sequential blocked schedule for the implemented model.

**Step 2: Isolate approximation language**

Reserve “approximate” for `independent_blocks=true`, which freezes the sweep-level residual, except where off-block weighted genotype crossproducts vanish.

**Step 3: Align Annotated BayesR wording**

Keep the existing distinction that annotated block BayesR uses the same marker-specific priors but a different transition schedule. Link to the newly titled BayesR3 strategy page and use its terminology consistently.

**Step 4: Search for stale approximation claims**

Run:

```bash
rg -n -i "accelerated approximation|heuristic approximation|exact sequential|independent-block" docs/src/manual
```

Expected: no supported sequential block path is called a heuristic or accelerated approximation; independent blocks remain clearly qualified.

**Step 5: Commit the exactness terminology**

```bash
git add docs/src/manual/block_bayesc.md docs/src/manual/annotated_bayesr.md
git commit -m "docs: separate sequential and independent block schedules"
```

### Task 5: Separate operation counts from the BayesR3 runtime fit

**Files:**
- Modify: `docs/src/manual/block_bayesc.md:221-299`

**Step 1: Retain the JWAS operation-count derivation**

Keep the derivation in terms of `N`, `P`, `s_i`, and requested chain length, labeling any ratio derived from it as an arithmetic operation-count comparison rather than measured speedup.

**Step 2: Retain the BayesR3 empirical evidence separately**

Describe `(N + b)/b` as a published empirical runtime fit from the BayesR3 implementation, datasets, libraries, hardware, and tested range.

**Step 3: Remove the pseudo-benchmark comparison**

Remove the shared coefficient table and claims such as:

```markdown
BayesR3 paper fit vs standard: `~446x` lower
the apparent `~2.0x` gap
```

Replace them with a direct statement that the two expressions must not be used as a software-to-software benchmark.

**Step 4: Review the performance section**

Run:

```bash
sed -n '235,330p' docs/src/manual/block_bayesc.md
```

Expected: the page contains useful operation-count and empirical-fit information without a direct numerical JWAS-versus-BayesR3 speed comparison.

**Step 5: Commit the performance clarification**

```bash
git add docs/src/manual/block_bayesc.md
git commit -m "docs: separate block cost model from BayesR3 timing fit"
```

### Task 6: Align navigation and cross-references

**Files:**
- Modify: `docs/make.jl:24`
- Modify: `docs/src/manual/workflow.md:23,317`
- Modify: `docs/src/manual/annotated_bayesc.md:150`
- Modify: `docs/src/manual/annotated_bayesr.md:354`
- Modify: `docs/src/manual/multitrait_annotated_bayesc.md:159`
- Modify: `docs/src/manual/bayesc_bayesr_comparison.md:152`
- Modify: `docs/src/manual/sem.md:221`

**Step 1: Rename the navigation label**

Change the Documenter navigation label from `Block BayesC` to `BayesR3 Block Strategy` while keeping the target `manual/block_bayesc.md`.

**Step 2: Update cross-reference labels**

Change visible link text from `Block BayesC` to `BayesR3 Block Strategy` where the target page is being cited for general block behavior. Preserve the existing relative link target.

**Step 3: Correct workflow scope**

Replace the claim that `fast_blocks` enables only the block BayesC path with wording that it enables the generalized BayesR3 strategy for the supported model paths.

**Step 4: Search for stale page labels**

Run:

```bash
rg -n "Block BayesC|block BayesC path|block_bayesc" docs/make.jl docs/src
```

Expected: remaining uses of “block BayesC” describe the actual BayesC sampler, while page links and navigation use “BayesR3 Block Strategy.”

**Step 5: Commit navigation alignment**

```bash
git add docs/make.jl docs/src/manual
git commit -m "docs: align BayesR3 strategy navigation"
```

### Task 7: Build, inspect, and record the documentation change

**Files:**
- Create: `docs/plans/2026-07-10-generalized-bayesr3-documentation-implementation.md`
- Verify: `docs/build/manual/block_bayesc/index.html`

**Step 1: Check Markdown and terminology**

Run:

```bash
git diff --check
rg -n -i "JWAS uses BayesC \(not BayesR\)|accelerated approximation|BayesR3 paper fit vs standard|apparent.*gap" docs/src docs/make.jl
```

Expected: `git diff --check` succeeds and the stale phrases are absent.

**Step 2: Build the documentation**

Run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected: Documenter completes successfully without missing-page or unresolved-link errors.

**Step 3: Inspect the rendered output**

Run:

```bash
rg -n "BayesR3 Block Strategy|user-defined|unequal|empirical runtime fit|independent_blocks" docs/build/manual/block_bayesc/index.html
```

Expected: the rendered page contains the new title and key distinctions.

**Step 4: Write the implementation record**

Create `docs/plans/2026-07-10-generalized-bayesr3-documentation-implementation.md` with:

- scope and terminology decisions;
- files changed;
- confirmation that sampler behavior and API were unchanged;
- verification commands and outcomes;
- any remaining limitations, especially contiguous-only explicit blocks.

**Step 5: Review the final patch**

Run:

```bash
git status --short
git diff --stat HEAD
git diff -- docs/make.jl docs/src/manual docs/plans/2026-07-10-generalized-bayesr3-documentation-implementation.md
```

Expected: only the intended documentation files and implementation record are included; unrelated user files remain untouched.

**Step 6: Commit the verified documentation**

```bash
git add docs/make.jl docs/src/manual \
  docs/plans/2026-07-10-generalized-bayesr3-documentation-implementation.md
git commit -m "docs: document generalized BayesR3 strategy"
```
