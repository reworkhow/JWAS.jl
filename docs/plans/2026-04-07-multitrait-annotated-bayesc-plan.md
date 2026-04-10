# Multi-Trait Annotated BayesC Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add dense 2-trait annotated BayesC to JWAS using a 3-step tree annotation prior that reconstructs a coherent marker-specific joint prior over `00`, `10`, `01`, and `11`.

**Architecture:** Keep the current dense multi-trait BayesC marker-effect and variance updates. Generalize the annotation container so 2-trait annotated BayesC becomes a 3-step annotation case that rebuilds `snp_pi[m, 4]` after each marker sweep, then feed those marker-specific state priors into the existing multi-trait BayesC samplers.

**Tech Stack:** Julia, JWAS dense MCMC path, Distributions.jl, DataFrames.jl, Documenter.jl

---

### Task 1: Add failing tests for dense 2-trait annotated BayesC setup

**Files:**
- Modify: `test/unit/test_annotated_bayesc.jl`
- Modify: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Write setup and validation tests**

Add tests that cover:

- dense 2-trait `method="BayesC"` accepts `annotations=...`
- 2-trait annotated BayesC still rejects:
  - `storage=:stream`
  - `fast_blocks=true`
  - `constraint=true`
- the initialized annotation container for 2 traits has:
  - `nsteps = 3`
  - `nclasses = 4`
  - `coefficients` with 3 columns
  - `snp_pi` with 4 columns
- the initial `snp_pi` rows reproduce the supplied starting joint prior

**Step 2: Write a failing dense run test**

Add a small dense 2-trait annotated BayesC run that currently fails because the production path still rejects annotations or does not route them into the multi-trait sampler.

**Step 3: Run tests to verify the failure**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl"); include("test/unit/test_multitrait_mcmc.jl")'
```

Expected:

- the new setup tests fail on the current single-trait-only annotation plumbing

### Task 2: Refactor annotation setup so 2-trait BayesC can build a 3-step annotation state

**Files:**
- Modify: `src/1.JWAS/src/markers/readgenotypes.jl`
- Modify: `src/1.JWAS/src/build_MME.jl`
- Modify: `src/1.JWAS/src/types.jl`
- Modify: `src/1.JWAS/src/input_data_validation.jl`

**Step 1: Delay multi-trait BayesC annotation specialization until trait count is known**

Keep annotation file parsing and QC in `readgenotypes.jl`, but stop forcing BayesC annotations into the single-trait scalar/vector `Pi` path before `ntraits` is available.

**Step 2: Build a 3-step annotation container for 2 traits**

Once `ntraits == 2` is known, create annotation state with:

- `nsteps = 3`
- `nclasses = 4`
- zeroed coefficients and linear predictors
- `snp_pi` initialized by repeating the supplied starting joint `Pi` row across markers

**Step 3: Preserve existing single-trait behavior**

Keep the single-trait annotated BayesC initialization unchanged as the 1-step case.

**Step 4: Tighten unsupported-path checks**

Reject multi-trait annotated BayesC when:

- `storage != :dense`
- `fast_blocks != false`
- `constraint == true`
- `ntraits != 2`

### Task 3: Implement the 3-step annotation Gibbs update for 2-trait BayesC

**Files:**
- Modify: `src/1.JWAS/src/MCMC/annotation_updates.jl`
- Test: `test/unit/test_annotated_bayesc.jl`

**Step 1: Add 2-trait tree indicator helpers**

Derive:

- `z1_j = 1(delta_j != 00)` on all markers
- `z2_j = 1(delta_j == 11)` on markers in `{10, 01, 11}`
- `z3_j = 1(delta_j == 10)` on markers in `{10, 01}`

**Step 2: Reuse the current binary probit machinery step-by-step**

Implement a 3-step annotation updater that:

- subsets markers for each step
- samples liabilities and coefficients for each step
- rebuilds the 4-state `snp_pi` row for each marker

**Step 3: Keep single-trait BayesC unchanged**

Do not regress the current 1-step BayesC annotation update.

### Task 4: Hook per-marker 4-state priors into the dense multi-trait BayesC samplers

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Add marker-specific prior lookup for sampler I**

Replace shared `log(BigPi[state])` lookups with marker-row-specific log priors from `annotations.snp_pi[j, :]` for the neighboring states being compared.

**Step 2: Add marker-specific prior lookup for sampler II**

Use the same `annotations.snp_pi[j, :]` row as the full-state prior when sampler II is active.

**Step 3: Update the MCMC order**

For supported 2-trait annotated BayesC runs, after the multi-trait marker sweep:

1. update the annotation state
2. refresh the marker-level joint prior
3. continue with the existing variance updates

### Task 5: Update multi-trait prior helpers and summaries

**Files:**
- Modify: `src/1.JWAS/src/markers/tools4genotypes.jl`
- Modify: `src/1.JWAS/src/output.jl`
- Modify: `src/1.JWAS/src/JWAS.jl`

**Step 1: Support marker-level 4-state startup priors where needed**

Make the multi-trait BayesC startup helpers accept the annotation-backed per-marker 4-state prior where prior-to-marker variance setup depends on inclusion probabilities.

**Step 2: Keep user-facing `Pi` summaries interpretable**

For annotated 2-trait BayesC, summarize `Pi` as the average marker-level prior over the four states while keeping the per-marker matrix internal to sampling.

**Step 3: Add step-aware annotation-coefficient output**

Write the three 2-trait annotation steps clearly in output files and on-screen summaries.

### Task 6: Update documentation

**Files:**
- Modify: `docs/src/manual/annotated_bayesc.md`
- Modify: `docs/src/manual/bayesc_bayesr_comparison.md`
- Modify: `docs/src/manual/workflow.md`

**Step 1: Document the 2-trait annotated BayesC prior**

Explain:

- the 3-step tree parameterization
- the induced 4-state joint prior
- the dense-only scope

**Step 2: Document startup and update order**

State explicitly that:

- the supplied joint `Pi` seeds the first sweep
- annotation coefficients start at zero
- annotation fitting begins after the first marker sweep

### Task 7: Verification

**Files:**
- None

**Step 1: Run targeted tests**

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

**Step 2: Run broader regression tests**

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_advanced_coverage.jl")'
julia --project=. --startup-file=no test/runtests.jl
```

**Step 3: Run docs build if docs were updated**

```bash
julia --project=docs --startup-file=no docs/make.jl
```
