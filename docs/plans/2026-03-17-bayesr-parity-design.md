# BayesR Parity Benchmark And Review Design

## Goal

Validate the new single-trait dense BayesR implementation against the R reference at the posterior-summary level, then use those results to drive a focused review of the BayesR code path.

## Scope

This design covers:

- a reproducible parity harness using one shared synthetic dataset
- comparison of JWAS BayesR and the R reference on the same exported files
- summary-level parity criteria, not draw-by-draw identity
- a structured implementation review tied to benchmark results

This design does not cover:

- `storage=:stream`
- `fast_blocks`
- annotations
- user-configurable `gamma` or starting `pi`
- CI-hardening of the benchmark

## Approach Options

### Option 1: One Shared Dataset, Two Runners

Generate one synthetic dataset once, export it, run both JWAS BayesR and the R reference on the same files, then compare posterior summaries.

Pros:

- makes preprocessing differences explicit
- produces a reusable artifact
- strongest parity signal

Cons:

- needs some file plumbing
- requires a small comparison workflow

Recommendation: use this option.

### Option 2: Recreate The Simulation Separately In Both Languages

Generate the same data independently in Julia and R and compare results.

Pros:

- fewer exported files

Cons:

- easy to drift on preprocessing or encoding
- parity failures become ambiguous

### Option 3: Compare Only Inside Julia

Port the R logic into Julia and compare the two Julia paths.

Pros:

- fastest to build

Cons:

- not a real language-to-language parity check
- can mask the same assumption on both sides

## Design

### 1. Parity Harness

Create a dedicated benchmark/parity harness under `benchmarks/` rather than reusing `benchmarks/jwas_full_benchmark.jl`.

The harness should contain:

- a Julia dataset export and JWAS runner
- an R runner for the reference BayesR script
- a comparison step that reads both summary outputs and reports differences

The first parity artifact should be one small single-trait dense synthetic dataset with explicit marker ordering, phenotype ordering, and hyperparameter settings.

### 2. Matching Rules

The parity comparison must align:

- genotype matrix after preprocessing
- phenotype vector
- `gamma = [0.0, 0.01, 0.1, 1.0]`
- starting `pi`
- chain length and burn-in
- marker ordering
- interpretation of the shared BayesR variance

The parity comparison should not require:

- identical RNG streams across Julia and R
- identical class assignments iteration by iteration
- machine-precision agreement

Success criterion: summary-level parity.

### 3. Comparison Modes

Run two modes:

- `estimatePi=false`
- `estimatePi=true`

The fixed-`pi` mode isolates the effect sampler and variance updates. The estimated-`pi` mode then checks the Dirichlet update path separately.

### 4. Summary Outputs

Compare these posterior summaries:

- posterior mean `pi` vector
- posterior mean shared variance `sigmaSq`
- posterior mean residual variance
- posterior mean nonzero marker frequency
- posterior mean marker effects

For marker effects, compare:

- full-vector correlation
- absolute differences on the largest-effect markers

### 5. Tolerances

Use practical summary-level tolerances rather than exact equality.

Initial tolerance style:

- `pi`: per-class absolute tolerance
- shared variance and residual variance: relative tolerance
- marker effects: high correlation plus top-marker absolute-difference checks
- nonzero frequency: aggregate and selected per-marker absolute tolerance

These tolerances should start as benchmark/report thresholds, not CI-gating unit-test thresholds.

### 6. Review Workflow

Review the BayesR implementation alongside the parity results in this order:

1. `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
2. `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
3. `src/1.JWAS/src/markers/tools4genotypes.jl`
4. `src/1.JWAS/src/variance_components.jl`
5. `src/1.JWAS/src/output.jl`

Review focus:

- conditional-update math
- zero-class handling
- log-probability normalization
- `yCorr` bookkeeping
- initialization and `pi` updates
- prior scale derivation and shared variance interpretation
- output semantics for BayesR

### 7. Reporting

Write a short benchmark note under `benchmarks/reports/` summarizing:

- what matched
- what differed
- whether differences look expected, numerical, or like implementation bugs
- follow-up actions if needed

## Deliverables

- one BayesR parity dataset export path
- one Julia parity runner
- one R parity runner
- one comparison script or comparison step
- one short parity report
- one structured review pass of the BayesR implementation using the benchmark results
