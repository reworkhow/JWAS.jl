# SEM Documentation Page Design

**Date:** 2026-02-26

## Goal

Create one dedicated SEM documentation page in the Manual that serves three audiences in sequence on a single page:
1. Beginner users (how to run it correctly)
2. Reference users (exact argument/output behavior)
3. Advanced users/contributors (implementation and math flow)

## Scope

### In scope
- Add one new Manual page: `docs/src/manual/sem.md`
- Add a corresponding Manual navigation entry in `docs/make.jl`
- Keep SEM coverage in `workflow.md` concise and add a pointer to the new SEM page to reduce duplication drift
- Include practical troubleshooting notes for common SEM setup/runtime errors

### Out of scope
- Changing SEM algorithms or runtime behavior
- Reworking all existing workflow material; only minimal cross-link adjustments
- Splitting SEM docs into multiple pages (this design intentionally keeps one layered page)

## Information Architecture (Single Page)

Page title: **SEM: Beginner to Advanced**

Planned section order:
1. Intro scope paragraph
2. Quick index (jump links)
3. Beginner Quick Start
4. Reference
5. Deep Technical Notes
6. Troubleshooting
7. Related Pages

Rationale: preserves a progressive learning path while still allowing fast access through jump links.

## Section Design

### 1) Beginner Quick Start

Purpose: get a user from zero to a correct SEM run.

Content:
- What SEM in JWAS does (plain-language summary)
- Prerequisites and constraints:
  - Multi-trait only
  - `causal_structure` lower-triangular requirement
  - Column index affects row index interpretation
  - SEM path behavior that enforces residual covariance constraint and disables missing phenotype imputation
- End-to-end runnable example aligned with current behavior
- Expected SEM-related outputs and where they appear
- Quick sanity checks after running

### 2) Reference

Purpose: fast lookup for exact API/runtime behavior.

Content:
- `runMCMC(...; causal_structure=...)` behavior table
- Validation rules and side effects table
- `causal_structure` semantics with small matrix examples
- Output file reference table (file, trigger, meaning)
- Compatibility/limitations table

### 3) Deep Technical Notes

Purpose: bridge user docs and contributor-level understanding.

Content:
- Notation used by implementation (`Y`, `Λ`, `Λy`)
- Construction flow and key functions
- MCMC call-site flow for sampling structural coefficients
- Marker-effect post-processing pipeline (direct/indirect/overall)
- Practical caveats/assumptions observed in current code
- Pointers to source files/functions

### 4) Troubleshooting

Purpose: rapid diagnosis for common setup errors.

Content:
- Lower-triangular structure errors
- Single-trait misuse
- Naming mismatch around genotype-term outputs
- Version note for the resolved issue #162 regression

### 5) Related Pages

- Workflow
- Public API
- Block BayesC / Memory pages as contextual references where relevant

## Writing Style and UX Decisions

- One page only, layered by user maturity
- Keep beginner section executable with copy-paste blocks
- Keep reference section table-heavy and concise
- Keep deep section implementation-accurate but still readable
- Prefer internal links to existing manual pages over duplicated long explanations

## Verification Plan

1. Build docs locally:
   - `julia --project=docs docs/make.jl`
2. Confirm new page appears in Manual navigation.
3. Confirm anchors/quick index links render and navigate correctly.
4. Spot-check references and code snippets for current API names.

## Risks and Mitigations

- Risk: page becomes too long.
  - Mitigation: strict section boundaries + quick index at top.
- Risk: beginner and deep sections drift over time.
  - Mitigation: keep shared facts centralized in reference tables and link from other docs.
- Risk: workflow duplicates SEM instructions.
  - Mitigation: minimize SEM duplication in workflow and point to SEM page.
