# Multi-Trait BayesC Sampler Selection Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add an explicit `multi_trait_sampler` option so multi-trait BayesC, including annotated 2-trait BayesC, can intentionally use Gibbs sampler I or II while preserving current `:auto` behavior by default.

**Architecture:** Store the sampler preference on `Genotypes`, expose it through `get_genotypes`, validate supported combinations during model checks, and let `MTBayesABC!` dispatch on `:auto`, `:I`, or `:II`. Keep the existing sampler kernels and only change the selection logic.

**Tech Stack:** Julia, JWAS marker/MCMC code paths, unit tests in `test/unit`

---

### Task 1: Add the public sampler-selection option

**Files:**
- Modify: `src/1.JWAS/src/types.jl`
- Modify: `src/1.JWAS/src/markers/readgenotypes.jl`
- Test: `test/unit/test_annotated_bayesc.jl`

**Step 1: Write the failing test**

Add a small construction test that creates annotated 2-trait BayesC genotypes with:

```julia
geno = get_genotypes(
    genofile,
    G;
    method="BayesC",
    quality_control=false,
    annotations=annotations,
    Pi=start_pi,
    multi_trait_sampler=:II,
)

@test geno.multi_trait_sampler == :II
```

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: FAIL because `get_genotypes` does not accept `multi_trait_sampler` yet or `Genotypes` has no such field.

**Step 3: Write minimal implementation**

- Add `multi_trait_sampler` to `Genotypes` in `src/1.JWAS/src/types.jl`.
- Initialize it to `:auto` in the inner constructor.
- Add keyword `multi_trait_sampler::Symbol = :auto` to `get_genotypes` in `src/1.JWAS/src/markers/readgenotypes.jl`.
- Validate `multi_trait_sampler in (:auto, :I, :II)` in `get_genotypes`.
- Store the selected symbol on the returned `Genotypes` object for both `:dense` and `:stream` construction paths.
- Update the `get_genotypes` docstring bullets to describe `:auto`, `:I`, and `:II`.

**Step 4: Run test to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: PASS for the new construction/field test.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/types.jl src/1.JWAS/src/markers/readgenotypes.jl test/unit/test_annotated_bayesc.jl
git commit -m "feat: add multi-trait sampler selection option"
```

### Task 2: Validate supported combinations explicitly

**Files:**
- Modify: `src/1.JWAS/src/input_data_validation.jl`
- Test: `test/unit/test_annotated_bayesc.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Write the failing tests**

Add validation coverage for:

```julia
@test_throws ErrorException build_model(... annotated geno with multi_trait_sampler=:bogus ...)
@test_throws ErrorException build_model(... annotated multi-trait BayesC with multi_trait_sampler=:II and fast_blocks=true ...)
```

and one allowed case:

```julia
model = build_model(... annotated 2-trait BayesC using geno with multi_trait_sampler=:II ...)
@test model.M[1].multi_trait_sampler == :II
```

**Step 2: Run tests to verify they fail**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: FAIL because validation does not yet understand the new option.

**Step 3: Write minimal implementation**

- In `src/1.JWAS/src/input_data_validation.jl`, add explicit checks that:
  - `multi_trait_sampler` is one of `:auto`, `:I`, `:II`
  - forcing `:II` does not relax existing annotated multi-trait BayesC restrictions
  - unsupported methods/paths still error clearly if `:I`/`:II` is requested outside the multi-trait BayesC context
- Keep `:auto` as the compatibility default.

**Step 4: Run tests to verify they pass**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/input_data_validation.jl test/unit/test_annotated_bayesc.jl test/unit/test_multitrait_mcmc.jl
git commit -m "test: validate multi-trait sampler selection"
```

### Task 3: Route `MTBayesABC!` through the requested sampler

**Files:**
- Modify: `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl`
- Test: `test/unit/test_multitrait_mcmc.jl`

**Step 1: Write the failing test**

Add a focused annotated 2-trait BayesC smoke test that sets:

```julia
multi_trait_sampler=:II
```

and verifies the run completes and still emits:

```julia
@test haskey(output, "annotation coefficients annotated_mt_geno")
@test haskey(output, "pi_annotated_mt_geno")
```

This should be a true runtime test, not just construction.

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: FAIL or still exercise sampler I because dispatch ignores the new override.

**Step 3: Write minimal implementation**

- Update `MTBayesABC!` in `src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl` so sampler selection becomes:
  - `:auto` -> existing rule based on `length(genotypes.π) == 2^nModels`
  - `:I` -> force `_MTBayesABC_samplerI!`
  - `:II` -> force `_MTBayesABC_samplerII!`
- Keep prior-source selection (`GlobalPiPrior` vs `MarkerSpecificPiPrior`) unchanged.
- Keep `state_labels = collect(keys(genotypes.π))` explicit for forced sampler II.
- Add a short comment noting that annotated 2-trait BayesC can now opt into sampler II explicitly while `:auto` preserves legacy behavior.

**Step 4: Run tests to verify it passes**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
```

Expected: PASS, including the new annotated `:II` smoke run.

**Step 5: Commit**

```bash
git add src/1.JWAS/src/markers/BayesianAlphabet/MTBayesABC.jl test/unit/test_multitrait_mcmc.jl test/unit/test_annotated_bayesc.jl
git commit -m "feat: support explicit multi-trait BayesC sampler selection"
```

### Task 4: Final verification and notes

**Files:**
- Modify: `docs/plans/2026-04-09-multitrait-bayesc-sampler-selection-plan.md`
- Create: `docs/plans/2026-04-09-multitrait-bayesc-sampler-selection-implementation.md`

**Step 1: Run focused verification**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'
julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'
```

Expected: PASS.

**Step 2: Run broader regression**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: PASS.

**Step 3: Write implementation note**

Record:
- final public API
- validation behavior
- default compatibility semantics for `:auto`
- the fact that annotated 2-trait BayesC can now intentionally use sampler II

in:

```text
docs/plans/2026-04-09-multitrait-bayesc-sampler-selection-implementation.md
```

**Step 4: Commit**

```bash
git add docs/plans/2026-04-09-multitrait-bayesc-sampler-selection-plan.md docs/plans/2026-04-09-multitrait-bayesc-sampler-selection-implementation.md
git commit -m "docs: record multi-trait BayesC sampler selection work"
```
