# Memory Guardrails for Large BayesC Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add fail-fast/warn/off memory feasibility guardrails for marker allocations in `runMCMC` so infeasible dense BayesC runs are blocked early with actionable diagnostics.

**Architecture:** Introduce a small marker-memory estimator utility in `tools4genotypes.jl`, plus a guard policy check function in the same module. Wire these into `runMCMC` before marker-heavy setup. Keep BayesC math and existing allocation paths unchanged.

**Tech Stack:** Julia, JWAS internals (`runMCMC`, marker tools), `Test` stdlib.

---

### Task 1: Add Failing Unit Tests for Memory Estimation and Guard Policy

**Files:**
- Create: `test/unit/test_memory_guardrails.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Create tests that expect the following API to exist and behave:
- `JWAS.estimate_marker_memory(...)`
- `JWAS.check_marker_memory_guard!(...)`

Test cases:
1. Non-block + unit weights estimates only `X` (+small vectors)
2. Non-block + non-unit weights includes extra `N*P`
3. Block mode includes `XRinvArray` and `XpRinvX` terms
4. Guard mode `:error` throws when estimate exceeds threshold
5. Guard mode `:warn` does not throw
6. Guard mode `:off` bypasses checks
7. Invalid mode / invalid ratio throw clear errors

**Step 2: Run test to verify it fails**

Run:
```bash
cd /Users/haocheng/Github/JWAS.jl
julia --project -e 'using Test, JWAS; include("test/unit/test_memory_guardrails.jl")'
```

Expected: FAIL due to missing functions.

### Task 2: Implement Minimal Estimator and Guard Policy Helpers

**Files:**
- Modify: `src/1.JWAS/src/markers/tools4genotypes.jl`

**Step 1: Write minimal implementation**

Add:
- `estimate_marker_memory(nObs, nMarkers; element_bytes, has_nonunit_weights, block_starts)` returning a named tuple with:
  - `bytes_X`
  - `bytes_xRinvArray`
  - `bytes_XRinvArray`
  - `bytes_XpRinvX`
  - `bytes_xpRinvx`
  - `bytes_total`
- `format_bytes_human(bytes)` helper for readable messages
- `check_marker_memory_guard!(; mode, ratio, estimated_bytes, total_memory_bytes, context_string)`

Use integer-safe arithmetic (`Int128`) for large dimensions.

**Step 2: Run focused test to verify it passes**

Run:
```bash
cd /Users/haocheng/Github/JWAS.jl
julia --project -e 'using Test, JWAS; include("test/unit/test_memory_guardrails.jl")'
```

Expected: PASS for estimator and guard policy tests.

### Task 3: Wire Guardrails Into `runMCMC`

**Files:**
- Modify: `src/1.JWAS/src/JWAS.jl`

**Step 1: Write failing integration-style test hook**

Extend `test/unit/test_memory_guardrails.jl` to call `runMCMC` with intentionally strict settings on tiny demo data and confirm:
- `memory_guard=:error` with tiny `memory_guard_ratio` throws before MCMC
- `memory_guard=:off` runs normally on same tiny data

**Step 2: Run test to verify it fails**

Run:
```bash
cd /Users/haocheng/Github/JWAS.jl
julia --project -e 'using Test, JWAS; include("test/unit/test_memory_guardrails.jl")'
```

Expected: FAIL because `runMCMC` does not yet accept these keywords.

**Step 3: Write minimal implementation**

In `runMCMC`:
- Add keyword args:
  - `memory_guard = :error`
  - `memory_guard_ratio::Float64 = 0.80`
- After block setup (`fast_blocks`) and before heavy marker setup:
  - Detect precision (`double_precision ? 8 : 4` bytes)
  - Detect non-unit weights from `mme.invweights` after `getMME` is initialized path (or equivalent available state)
  - Build block starts when applicable
  - Call estimator
  - Call guard checker

Message should include estimated total and major components.

**Step 4: Run test to verify it passes**

Run:
```bash
cd /Users/haocheng/Github/JWAS.jl
julia --project -e 'using Test, JWAS; include("test/unit/test_memory_guardrails.jl")'
```

Expected: PASS.

### Task 4: Register Test and Run Project Test Subset

**Files:**
- Modify: `test/runtests.jl`

**Step 1: Add unit test include**

Include `test/unit/test_memory_guardrails.jl` in the unit test section.

**Step 2: Run lightweight verification**

Run:
```bash
cd /Users/haocheng/Github/JWAS.jl
julia --project -e 'using Test, JWAS; include("test/unit/test_memory_guardrails.jl")'
```

Expected: PASS.

**Step 3: Run standard suite command used in repo (if feasible)**

Run:
```bash
cd /Users/haocheng/Github/JWAS.jl
julia --project test/runtests.jl
```

Expected: PASS (or report exact failures if unrelated/time-bound).

### Task 5: Commit

**Files:**
- Add/modify all files above

**Step 1: Commit changes**

```bash
cd /Users/haocheng/Github/JWAS.jl
git add src/1.JWAS/src/JWAS.jl src/1.JWAS/src/markers/tools4genotypes.jl test/unit/test_memory_guardrails.jl test/runtests.jl docs/plans/2026-02-23-memory-guardrails-implementation.md
git commit -m "feat: add memory feasibility guardrails for large marker runs"
```
