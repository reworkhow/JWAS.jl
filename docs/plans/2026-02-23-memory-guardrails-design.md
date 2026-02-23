# Path 1 Design: Memory Feasibility Guardrails for Original BayesC

## Goal
Prevent out-of-memory failures for large marker analyses by estimating marker-path memory before heavy allocations and failing fast (with explicit override) when the estimate exceeds available system memory.

## Scope
- In scope: `runMCMC` guardrail behavior for marker memory feasibility, including `fast_blocks=false` and `fast_blocks!=false` cases.
- Out of scope: changing BayesC sampling math, changing genotype representation, out-of-core backends.

## Design Summary
Implement a pre-allocation memory estimator and policy check before MCMC starts marker-matrix setup.

1. Add an estimator function in `src/1.JWAS/src/markers/tools4genotypes.jl` that computes:
   - Base dense genotype memory (`X`: `N*P*bytes_per_value`)
   - Extra weighted copy for non-unit residual weights (`xRinvArray`)
   - Block extras when `fast_blocks != false`:
     - `XRinvArray` (`N*P` values)
     - `XpRinvX` (`sum_i s_i^2` values)
2. Add policy controls in `runMCMC(...)`:
   - `memory_guard=:error` (default)
   - `memory_guard_ratio=0.80` (fraction of total system RAM allowed)
3. Before heavy marker preparation, evaluate estimated memory against threshold:
   - If below threshold: continue silently
   - If above threshold:
     - `:error`: throw with detailed, actionable message
     - `:warn`: print warning and continue
     - `:off`: skip check
4. Keep current BayesC behavior unchanged except preflight guardrails.

## Data Flow
1. User calls `runMCMC(...)`.
2. `runMCMC` normal preprocessing runs.
3. Guardrail computes memory estimate from current marker dimensions/precision and block setting.
4. Guardrail policy decides fail/warn/continue.
5. Existing allocation and MCMC flow remain unchanged.

## Error Handling
- Invalid `memory_guard` value => immediate `error` with accepted values.
- Invalid `memory_guard_ratio` (`<=0` or `>1`) => immediate `error`.
- Estimation overflow risk is handled using integer-safe arithmetic with `Int128` during estimation and conversion back to display strings.

## Testing Strategy (TDD)
- Add estimator-focused tests for:
  - non-block unit weights
  - non-block non-unit weights
  - block mode extras
- Add guard policy tests for:
  - `:error` mode throws when threshold exceeded
  - `:warn` mode does not throw
  - `:off` mode bypasses checks
- Keep tests small and deterministic by using synthetic dimensions and forced low thresholds.

## Success Criteria
- For infeasible large dense configurations (e.g., `N=500k, P=2M`), users get immediate actionable failure/warning before MCMC enters heavy marker allocations.
- Existing feasible workflows remain behaviorally unchanged when estimates are below threshold.
