# Fast-Block Multi-Trait BayesC Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add block-form sampler I and sampler II support for multi-trait BayesC so plain and annotated multi-trait BayesC both work with `fast_blocks`, while making `multi_trait_sampler=:I` the default and retaining explicit `:II` and `:auto`.

**Architecture:** Keep the dense multi-trait BayesC path unchanged and refactor the block path to mirror it. The new block wrapper will choose a prior source and sampler mode, then dispatch to block translations of the current dense sampler-I and sampler-II kernels. Validation and tests will be updated so annotated multi-trait BayesC can use `fast_blocks` through the new block path.

**Tech Stack:** Julia, JWAS multi-trait BayesC samplers, existing unit-test suite under `test/unit/`.

---

### Task 1: Change the default multi-trait sampler to `:I`

**Files:**
- Modify: `src/1.JWAS/src/markers/readgenotypes.jl`
- Modify: `src/1.JWAS/src/types.jl`
- Modify: `docs/src/manual/multitrait_annotated_bayesc.md`
- Test: `test/unit/test_annotated_bayesc.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Write the failing tests**

Add assertions that:
- multi-trait BayesC created without an explicit `multi_trait_sampler` now stores `:I`
- explicit `:auto` still stores `:auto`

**Step 2: Run the focused tests to verify failure**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: tests fail because the current default is still `:auto`.

**Step 3: Make the minimal implementation change**

Update the default value and any user-facing documentation strings that still describe `:auto` as the default.

**Step 4: Re-run the focused tests**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: the new default assertions pass.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/readgenotypes.jl src/1.JWAS/src/types.jl docs/src/manual/multitrait_annotated_bayesc.md test/unit/test_annotated_bayesc.jl test/unit/test_multitrait_mcmc.jl
git commit -m "refactor: default multitrait BayesC to sampler I"
```

### Task 2: Refactor the block wrapper to mirror the dense wrapper

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Write the failing tests**

Add direct low-level tests that the block wrapper:
- uses `GlobalPiPrior` for plain multi-trait BayesC
- uses `MarkerSpecificPiPrior` for annotated 2-trait BayesC
- respects explicit `multi_trait_sampler=:II`

Keep the tests narrow by probing dispatch behavior rather than full MCMC first.

**Step 2: Run the focused tests to verify failure**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: failure because the block wrapper currently always uses the old sampler-I-only path.

**Step 3: Write the wrapper refactor**

In `MTBayesABC.jl`:
- keep `mt_bayesc_sampler_mode(...)`
- make `MTBayesABC_block!(genotypes, ...)` choose the prior source like the dense wrapper
- make it choose sampler mode and dispatch to block sampler I or II

Do not change the dense wrapper.

**Step 4: Re-run the focused tests**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: new dispatch tests pass, but block sampler II is still missing.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl test/unit/test_multitrait_mcmc.jl
git commit -m "refactor: route block multitrait BayesC by sampler mode"
```

### Task 3: Implement block sampler I using the dense sampler-I logic

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`
- Test: `test/unit/test_annotated_bayesc.jl`

**Step 1: Write the failing tests**

Add focused block-run tests for:
- plain multi-trait BayesC with `fast_blocks=true`
- annotated multi-trait BayesC with `fast_blocks=true`

Use short chains and assert successful execution plus expected output shape.

**Step 2: Run the focused tests to verify failure**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: annotated block run still fails under validation or old block behavior.

**Step 3: Implement `_MTBayesABC_block_samplerI!`**

Translate the current dense sampler-I algebra into the block-local residual form:
- block-local `XpRinvycorr`
- `XpRinvX` updates
- same `log_prior(...)`
- same one-indicator-at-a-time state updates

Keep the old block body only long enough to verify the translation, then remove it.

**Step 4: Re-run the focused tests**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: plain and annotated block sampler-I runs pass.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl test/unit/test_multitrait_mcmc.jl test/unit/test_annotated_bayesc.jl
git commit -m "feat: add block sampler I for multitrait BayesC"
```

### Task 4: Implement block sampler II using the dense sampler-II logic

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`
- Test: `test/unit/test_annotated_bayesc.jl`

**Step 1: Write the failing tests**

Add short block-run tests for:
- plain multi-trait BayesC with `fast_blocks=true, multi_trait_sampler=:II`
- annotated multi-trait BayesC with `fast_blocks=true, multi_trait_sampler=:II`

**Step 2: Run the focused tests to verify failure**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: explicit block sampler-II runs fail or route incorrectly.

**Step 3: Implement `_MTBayesABC_block_samplerII!`**

Translate dense sampler-II logic into block-local algebra:
- precompute block-local state support
- use `log_prior(...)`
- use the same joint-state draw logic
- update `XpRinvycorr` after each sampled state

**Step 4: Re-run the focused tests**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: explicit `:II` block runs pass.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl test/unit/test_multitrait_mcmc.jl test/unit/test_annotated_bayesc.jl
git commit -m "feat: add block sampler II for multitrait BayesC"
```

### Task 5: Lift annotated multi-trait fast-block validation restrictions

**Files:**
- Modify: `src/1.JWAS/src/input_data_validation.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Test: `test/unit/test_annotated_bayesc.jl`

**Step 1: Write the failing tests**

Replace or update the current tests that expect annotated multi-trait BayesC to reject `fast_blocks=true`.

**Step 2: Run the focused test to verify failure**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: the old validation rejects the new supported path.

**Step 3: Make the validation changes**

Remove only the `fast_blocks=false` restriction for annotated multi-trait BayesC. Keep the other current restrictions:
- exactly 2 traits
- `storage=:dense`
- `constraint=false`

Adjust any runtime guard that still rejects annotated multi-trait BayesC in block mode.

**Step 4: Re-run the focused test**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: the supported annotated block cases now run.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/input_data_validation.jl src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl test/unit/test_annotated_bayesc.jl
git commit -m "feat: allow fast-block annotated multitrait BayesC"
```

### Task 6: Add broader regression coverage and benchmark smoke coverage

**Files:**
- Modify: `test/unit/test_multitrait_mcmc.jl`
- Modify: `test/unit/test_misc_coverage.jl`
- Optional Modify: `benchmarks/simulated_annotations_multitrait_comparison.jl`

**Step 1: Add the regression cases**

Cover:
- default multi-trait sampler now `:I`
- explicit `:auto` still available
- block plain MT BayesC `:auto`
- block annotated MT BayesC `:auto`
- explicit block `:II`

If benchmark case enumeration depends on the old default assumption, update it explicitly rather than relying on defaults.

**Step 2: Run the focused tests**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: all new regression coverage passes.

**Step 3: Commit**

```bash
git add test/unit/test_multitrait_mcmc.jl test/unit/test_misc_coverage.jl benchmarks/simulated_annotations_multitrait_comparison.jl
git commit -m "test: cover block multitrait BayesC sampler modes"
```

### Task 7: Verify the whole affected surface

**Files:**
- Verify only

**Step 1: Run focused tests**

Run:
```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected: PASS

**Step 2: Run the full suite**

Run:
```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: PASS

**Step 3: If docs changed, rebuild docs**

Run:
```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected: PASS

**Step 4: Create the implementation record**

Create:
- `docs/plans/2026-04-12-mt-bayesc-fast-blocks-implementation.md`

Include:
- files changed
- tests run
- remaining limitations

**Step 5: Final commit**

```bash
git add docs/plans/2026-04-12-mt-bayesc-fast-blocks-implementation.md
git commit -m "docs: record multitrait BayesC fast-block implementation"
```
