# BayesR Fast-Blocks Benchmark Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Keep production BayesR untouched while separating production-facing fast-block benchmarks from exploratory local-loop diagnostics.

**Architecture:** Trim the main BayesR fast-block benchmark to production comparisons only, move exploratory modes into a dedicated debug script, and update tests to reflect the cleaner split.

**Tech Stack:** Julia, CSV.jl, DataFrames.jl, JWAS test suite

---

### Task 1: Write the failing benchmark test

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`

**Step 1: Write the failing test**

Change the long-chain benchmark subprocess test so it expects only two methods in
`schedule_runs.csv`:

- `dense`
- `burnin_gated`

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:

- the long-chain benchmark mode test fails because `single_rep` is still present

### Task 2: Refactor benchmark scripts

**Files:**
- Modify: `benchmarks/bayesr_fast_blocks_parity.jl`
- Create: `benchmarks/debug/bayesr_fast_blocks_debug.jl`

**Step 1: Trim the main script**

Keep:

- production-facing helpers
- `fixed_hyperparameters`
- `dense_multiseed`
- `long_chain_schedule_comparison`

Remove from the main script:

- benchmark-local trace helpers
- benchmark-local single-rep helpers
- debug-only modes

**Step 2: Create debug script**

Move the exploratory local-loop machinery into:

- `benchmarks/debug/bayesr_fast_blocks_debug.jl`

It should own:

- local BayesR state setup
- local trace mode
- local single-rep mode

**Step 3: Re-run benchmark test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected:

- benchmark test passes

### Task 3: Verify full suite and document outcome

**Files:**
- Create: `docs/plans/2026-03-20-bayesr-fast-blocks-benchmark-cleanup-implementation.md`

**Step 1: Run full verification**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected:

- full suite passes

**Step 2: Write implementation note**

Summarize:

- what stayed in the production-facing benchmark
- what moved to debug
- what tests now cover

**Step 3: Commit**

Commit the cleanup as a benchmark refactor after verification.
