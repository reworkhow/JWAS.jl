# BayesR Fast Blocks Long-Chain Schedule Comparison Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a reproducible long-chain BayesR benchmark that compares dense BayesR, default-blocks single-rep BayesR, and the current production burnin-gated fast-block schedule.

**Architecture:** Extend the existing BayesR fast-block benchmark harness with one new multiseed comparison mode. Reuse the current local single-rep helper and the production `runMCMC` path for dense and burnin-gated runs.

**Tech Stack:** Julia, JWAS benchmark helpers, CSV, DataFrames, Test

---

### Task 1: Add a failing subprocess test for the new benchmark mode

**Files:**
- Modify: `test/unit/test_bayesr_parity.jl`

**Step 1: Write the failing test**

Add a subprocess test for:

```bash
JWAS_BAYESR_BLOCK_MODE=long_chain_schedule_comparison \
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl <tmpdir>
```

with small test values for chain length, burnin, seeds, and dataset size.

Assert the mode writes:

- `schedule_runs.csv`
- `schedule_pairwise_summary.csv`

**Step 2: Run it and verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected: fail because the mode is unsupported.

### Task 2: Implement the long-chain schedule comparison mode

**Files:**
- Modify: `benchmarks/bayesr_fast_blocks_parity.jl`

**Step 1: Add helper(s) for multiseed method summaries**

Collect per-seed summaries for:

- dense production BayesR
- local default-blocks single-rep BayesR
- production burnin-gated BayesR fast blocks

**Step 2: Add `long_chain_schedule_comparison` mode**

Write:

- `schedule_runs.csv`
- `schedule_pairwise_summary.csv`

### Task 3: Run targeted verification

**Step 1: Run the targeted benchmark tests**

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr_parity.jl")'
```

Expected: pass.

### Task 4: Run the actual long-chain benchmark

**Step 1: Execute the benchmark**

```bash
JWAS_BAYESR_BLOCK_MODE=long_chain_schedule_comparison \
JWAS_BAYESR_BLOCK_CHAIN_LENGTH=10000 \
JWAS_BAYESR_BLOCK_BURNIN=2000 \
JWAS_BAYESR_BLOCK_SEEDS=2026,2027,2028,2029,2030 \
julia --project=. --startup-file=no benchmarks/bayesr_fast_blocks_parity.jl <outdir>
```

**Step 2: Inspect the summary outputs**

Determine whether:

- `single_rep` remains close to `dense`
- `burnin_gated` remains meaningfully separated from `dense`

### Task 5: Record implementation and findings

**Files:**
- Create: `docs/plans/2026-03-20-bayesr-fast-blocks-long-chain-schedule-comparison-implementation.md`

**Step 1: Save**

Document:

- files changed
- verification commands
- benchmark outputs
- interpretation of the long-chain result
