# BayesR Fast Blocks Default-Blocks Single-Rep Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a benchmark-only BayesR diagnostic that uses the default block partition but forces one within-block sweep per outer iteration, then compare it against dense BayesR with free `pi` and free `sigmaSq`.

**Architecture:** Extend the existing `benchmarks/bayesr_fast_blocks_parity.jl` script with one new mode and a local Gibbs driver that reuses production BayesR update math. Keep the change out of the production sampler path.

**Tech Stack:** Julia, JWAS benchmark helpers, CSV, DataFrames, Test

---

### Task 1: Add a failing benchmark-mode test

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`

**Step 1: Write the failing test**

Add a test that launches:

```bash
JWAS_BAYESR_BLOCK_MODE=default_blocks_single_rep julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl <tmpdir>
```

and expects the benchmark outputs:

- `comparison_scalar_metrics.csv`
- `comparison_pi.csv`
- `comparison_marker_effects_top.csv`
- `runtime.csv`

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected: failure because the benchmark mode is unsupported.

### Task 2: Implement the local single-rep diagnostic mode

**Files:**
- Modify: `benchmarks/bayesr_fast_blocks_parity.jl`

**Step 1: Add a local free-hyperparameter Gibbs helper**

Implement a helper that:

- initializes a local BayesR state
- runs the outer Gibbs loop
- reuses production helper math for `pi`, `sigmaSq`, and residual variance updates
- uses:
  - `BayesR!` for dense
  - `BayesR_block!(..., iter, typemax(Int))` for single-rep default blocks

**Step 2: Add `default_blocks_single_rep` mode**

Make the new mode:

- use default block partition
- keep outer `chain_length` and `burnin` equal to dense
- compare dense local summaries against single-rep default-block summaries
- write the standard comparison CSVs

### Task 3: Run targeted verification

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`
- Modify: `benchmarks/bayesr_fast_blocks_parity.jl`

**Step 1: Run the targeted test file**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected: pass.

**Step 2: Run the new diagnostic benchmark**

Run:

```bash
JWAS_BAYESR_BLOCK_MODE=default_blocks_single_rep julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl /tmp/bayesr_fast_blocks_default_single_rep
```

Expected: benchmark completes and writes the comparison CSVs.

### Task 4: Record implementation

**Files:**
- Create: `docs/plans/2026-03-19-bayesr-fast-blocks-default-blocks-single-rep-implementation.md`

**Step 1: Save what was implemented**

Document:

- files changed
- verification commands
- diagnostic results
- interpretation of whether default blocks with one sweep match dense BayesR
