# Annotated BayesR Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add single-trait dense annotated BayesR to JWAS as an individual-level analogue of `sbayesrc.R`, with annotation-driven 4-class mixture priors and step-specific annotation-coefficient output.

**Architecture:** Generalize the existing annotation framework so Annotated BayesC remains the 1-step case and annotated BayesR becomes the 3-step case. Keep the dense BayesR marker sampler intact except for replacing the shared class prior with per-SNP class probabilities derived from the current annotation coefficients.

**Tech Stack:** Julia, JWAS dense MCMC path, Distributions.jl, DataFrames.jl, Documenter.jl

---

### Task 1: Add failing annotated BayesR API tests

**Files:**
- Create: `test/unit/test_annotated_bayesr.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing tests**

Add a new test file that covers:

- `get_genotypes(...; method="BayesR", annotations=annotations)` is accepted for dense single-trait use
- annotated BayesR still rejects:
  - `fast_blocks=true`
  - `storage=:stream`
  - multi-trait
  - `RRM`
- the default annotation intercept conversion reproduces:
  - `p1 = 0.05`
  - `p2 = 0.4`
  - `p3 = 0.25`
- the initialized BayesR annotation container has:
  - `coefficients` size `(p, 3)`
  - `snp_pi` size `(m, 4)`

Include the new file from `test/runtests.jl`.

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected:

- failure because BayesR annotations are still rejected and the step-aware annotation container does not exist yet

**Step 3: Commit**

Do not commit yet. This task establishes the failing test target for the next implementation tasks.

### Task 2: Generalize the annotation container and input plumbing

**Files:**
- Modify: `src/1.JWAS/src/types.jl`
- Modify: `src/1.JWAS/src/markers/readgenotypes.jl`
- Modify: `src/1.JWAS/src/input_data_validation.jl`
- Test: `test/unit/test_annotated_bayesr.jl`
- Regression: `test/unit/test_annotated_bayesc.jl`

**Step 1: Implement a step-aware `MarkerAnnotations`**

Refactor `MarkerAnnotations` so it can represent either:

- `nsteps = 1` for Annotated BayesC
- `nsteps = 3` for annotated BayesR

Concrete changes:

- make `coefficients`, `mean_coefficients`, and `mean_coefficients2` matrices
- make `variance` a vector of step-specific variances
- make `liability`, `mu`, `lower_bound`, and `upper_bound` matrices with one column per step
- add `snp_pi` to store the current per-SNP class probabilities
- keep `design_matrix` and `lhs`

**Step 2: Implement annotated BayesR initialization**

In `readgenotypes.jl`:

- allow `annotations` for `method="BayesR"`
- build `MarkerAnnotations(design_matrix; nsteps=3, nclasses=4)`
- initialize the 3 intercepts from the default BayesR prior:
  - `Phi^{-1}(0.05)`
  - `Phi^{-1}(0.40)`
  - `Phi^{-1}(0.25)`
- initialize all non-intercept coefficients to zero
- initialize `snp_pi` so every SNP starts with:
  - `[0.95, 0.03, 0.015, 0.005]`

**Step 3: Preserve BayesC behavior**

Keep Annotated BayesC on the same container with `nsteps = 1` and ensure:

- existing annotation validation stays intact
- existing Annotated BayesC API and output behavior do not regress

**Step 4: Re-run targeted annotation tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl"); include("test/unit/test_annotated_bayesr.jl")'
```

Expected:

- Annotated BayesC tests pass
- BayesR annotation API tests progress further, but BayesR MCMC/update tests still fail because the sampler is not integrated yet

### Task 3: Generalize the annotation Gibbs helpers

**Files:**
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Test: `test/unit/test_annotated_bayesr.jl`
- Regression: `test/unit/test_annotated_bayesc.jl`

**Step 1: Refactor single-step annotation helpers into step-aware helpers**

Replace the current single-step BayesC-only logic with step-aware helpers that can:

- initialize the annotation state for `nsteps = 1` or `3`
- compute truncated-normal bounds per step
- sample liabilities per step
- sample annotation coefficients per step
- sample one annotation-effect variance per step
- rebuild `snp_pi`

**Step 2: Implement the BayesR step-up indicator logic**

For annotated BayesR, derive:

- `z1_j = 1(delta_j > 1)`
- `z2_j = 1(delta_j > 2)`
- `z3_j = 1(delta_j > 3)`

Then fit:

- step 1 on all SNPs
- step 2 on SNPs with `z1 = 1`
- step 3 on SNPs with `z2 = 1`

**Step 3: Keep BayesC on the same code path**

For Annotated BayesC:

- use the same helper framework with `nsteps = 1`
- keep the current binary-inclusion logic

**Step 4: Re-run targeted annotation tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl"); include("test/unit/test_annotated_bayesr.jl")'
```

Expected:

- helper-level initialization tests pass
- annotated BayesR end-to-end sampler tests still fail until the BayesR dense sweep consumes `snp_pi`

### Task 4: Integrate annotated priors into the dense BayesR sampler

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Test: `test/unit/test_annotated_bayesr.jl`
- Regression: `test/unit/test_bayesr.jl`

**Step 1: Write the failing dense annotated BayesR run test**

In `test/unit/test_annotated_bayesr.jl`, add a small dense single-trait run that:

- loads BayesR genotypes with `annotations=...`
- runs `runMCMC(...)`
- checks that the run completes
- checks that `Model_Frequency` exists in the marker-effects output

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected:

- failure because BayesR still uses the shared `pi` vector rather than per-SNP class probabilities

**Step 3: Implement the dense BayesR prior hookup**

Update the dense BayesR sweep so that:

- if `Mi.annotations === false`, it keeps using the shared `Mi.π`
- if `Mi.annotations !== false`, it uses row `j` of `ann.snp_pi` as the class prior for SNP `j`

Update the single-trait MCMC loop so annotated BayesR does:

1. marker sweep with current `ann.snp_pi`
2. annotation update
3. `sigmaSq` update
4. residual-variance update

**Step 4: Re-run targeted BayesR tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl"); include("test/unit/test_annotated_bayesr.jl")'
```

Expected:

- BayesR regression tests pass
- annotated BayesR dense run test passes

### Task 5: Add output and manual documentation

**Files:**
- Modify: `src/1.JWAS/src/output.jl`
- Create: `docs/src/manual/annotated_bayesr.md`
- Modify: `docs/make.jl`
- Modify: `docs/src/index.md`
- Modify: `docs/src/manual/workflow.md`
- Test: `test/unit/test_annotated_bayesr.jl`

**Step 1: Write the failing output-shape test**

Add a test that the annotated BayesR output contains:

- `output["annotation coefficients <Mi.name>"]`

with columns:

- `Annotation`
- `Step`
- `Estimate`
- `SD`

and that the `Step` column contains all three BayesR step labels.

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected:

- failure because output is still using the existing Annotated BayesC 3-column format

**Step 3: Implement output and docs**

In `output.jl`:

- keep existing Annotated BayesC output unchanged
- add a long-table branch for annotated BayesR with one row per `(annotation, step)`

In docs:

- add `manual/annotated_bayesr.md`
- add it to `docs/make.jl`
- add it to `docs/src/index.md`
- add a short link from `docs/src/manual/workflow.md`

The manual page should explain:

- the 3 conditional probabilities
- the 4-class reconstruction
- initialization from the default BayesR prior
- the Gibbs sampler order
- current restrictions
- a minimal usage example

**Step 4: Build the docs**

Run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected:

- docs build passes

### Task 6: Full verification and feature record

**Files:**
- Create: `docs/plans/2026-03-27-annotated-bayesr-implementation.md`

**Step 1: Run the full test suite**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected:

- full suite passes

**Step 2: Write the implementation note**

Summarize:

- what changed in the annotation container
- how annotated BayesR maps to `sbayesrc.R`
- what support is included in v1
- what remains deferred
- which tests and docs were run

**Step 3: Commit**

Commit the feature in small reviewable chunks. Recommended sequence:

1. annotation container and validation
2. sampler integration
3. output and docs
4. final implementation note
