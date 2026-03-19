# Prototype-to-Production Benchmarking Guidelines

This page explains how benchmarking should be done when a demo script, prototype sampler, or reference implementation is being converted into production JWAS package code.

This page is a process guide. It is not a benchmark report.

- Use this page to decide **how** to benchmark a new method or backend.
- Save concrete benchmark runs and conclusions under `benchmarks/reports/`.

## When To Use This Page

Use this guideline when any of the following is true:

- a new method is being added to JWAS from demo or prototype code
- a reference implementation exists in another language and JWAS needs parity checks
- a new backend is being added, such as streaming, block updates, annotations, or a new data representation
- an existing method is being refactored and you need to show behavior did not change
- a short benchmark result looks suspicious and you need to decide whether it is a bug or Monte Carlo noise

## Core Principles

### 1. Benchmark the production path first

The final claim must be based on the real JWAS path that users will call, not only on a helper script or local prototype.

For JWAS this usually means benchmarking:

- `get_genotypes(...)`
- `build_model(...)`
- `runMCMC(...)`

Local replay code and custom benchmark loops are useful for diagnosis, but they do not replace benchmarking the production path.

### 2. Separate deterministic correctness from stochastic agreement

A stochastic mismatch is not automatically a bug.

Benchmarking should be staged:

1. deterministic setup checks
2. one-step update checks
3. short-chain debugging
4. long-chain or multiseed summary comparison

If these stages are mixed together, it becomes difficult to tell whether the problem is:

- wrong formulas
- wrong preprocessing
- different initialization
- different random-number streams
- or ordinary Monte Carlo variation

### 3. Hold the benchmark inputs constant before blaming the sampler

Before comparing two implementations, align:

- genotype matrix after preprocessing
- phenotype vector
- marker ordering
- priors
- starting values
- chain length
- burn-in
- what is estimated and what is fixed

If any of those differ, parity conclusions are weak.

### 4. Do not trust a single short MCMC run

For Bayesian methods, one short chain can easily create a false alarm.

Short chains are useful for:

- fast debugging
- smoke tests
- identifying gross mistakes

They are not sufficient for final parity conclusions. Final parity should use:

- longer chains
- multiple seeds
- or both

### 5. Save enough artifacts that someone else can rerun the benchmark

A benchmark is not complete if it only exists as terminal output.

Every serious benchmark should save:

- dataset or dataset-generation config
- chain settings
- seeds
- scalar summaries
- method-specific summaries such as `pi`
- comparison tables
- a short markdown conclusion

## Benchmark Types Required For JWAS

When moving from prototype to production, not every project needs the same amount of benchmarking. But the following categories are the default checklist.

### A. Smoke benchmark

Purpose:

- prove the new path runs end to end
- prove outputs exist
- catch obvious validation or integration mistakes

Typical scope:

- small dataset
- short chain
- one seed

Questions answered:

- does the method run?
- does the package write expected outputs?
- are unsupported combinations rejected cleanly?

### B. Deterministic setup benchmark

Purpose:

- prove preprocessing and initialization match the intended design

Check items:

- same marker count
- same ordering
- same priors
- same initial variances
- same initial class labels or inclusion indicators when relevant
- same fixed or estimated settings for `pi`, variance, scale, and annotations

This stage prevents wasted debugging on mismatched setup.

### C. One-step controlled replay benchmark

Purpose:

- verify the actual update formulas
- find the first mismatch location precisely

This is the preferred tool when a reference implementation exists.

For a one-step replay:

- export a shared initial state
- export a shared list of random draws
- force both implementations to consume the same draws
- compare every substep

Typical quantities to compare:

- conditional means
- conditional variances
- class probabilities
- chosen class or inclusion state
- sampled effect
- sufficient statistics such as `ssq` and `nnz`
- updated variance draw

This benchmark answers:

- are the formulas the same?
- are the random draws mapped into the formulas the same way?

It does **not** answer whether long-run posterior summaries match.

### D. Production-vs-local benchmark

Purpose:

- check whether the discrepancy lives inside JWAS itself before bringing in another language

Typical comparison:

- production `runMCMC(...)`
- benchmark-local or debug-local Julia implementation of the same algorithm

If those two disagree materially, the issue is inside the JWAS stack.

If those two agree closely, the remaining question is usually cross-language behavior or benchmark protocol, not production integration.

### E. Long-chain parity benchmark

Purpose:

- compare posterior summaries at a scale where Monte Carlo noise is smaller

Typical outputs:

- posterior mean marker variance scale
- posterior mean residual variance
- posterior sparsity or nonzero frequency
- posterior `pi`
- marker-effect correlation or top-effect comparisons

This is the main parity benchmark for a stochastic method.

### F. Multiseed parity benchmark

Purpose:

- distinguish implementation mismatch from ordinary between-chain variability

If one short chain shows a large difference, the next question should be:

- is this difference larger than ordinary chain-to-chain variation?

Multiseed parity should save:

- per-seed summary table
- aggregate mean and max differences
- between-seed scale summaries for each implementation

For BayesR, this turned out to be essential. A short single chain suggested a `sigmaSq` mismatch that largely disappeared when longer multiseed runs were used.

### G. Performance and scale benchmarks

Purpose:

- measure runtime and memory cost separately from correctness

Do not mix correctness and scaling claims into one number.

Performance benchmarks should record:

- hardware
- Julia version
- storage mode
- precision mode
- block settings
- dataset dimensions
- chain settings
- timed region

Scale benchmarks should answer questions like:

- does runtime scale linearly or superlinearly with markers?
- how much memory does dense mode require?
- when should streaming or blocks be preferred?

## Recommended Workflow

The following is the default JWAS workflow for converting prototype code into package code.

### Step 1: Define the production target

Write down exactly what is being productionized.

Examples:

- single-trait dense BayesR
- streaming BayesC with unit weights
- annotated BayesC with marker-specific inclusion priors

Also write down what is **not** in scope yet. This keeps the benchmark honest.

### Step 2: Build a shared benchmark dataset

Create one dataset generator or one exported dataset that all implementations will use.

The benchmark dataset should record:

- data dimensions
- preprocessing
- marker order
- seeds for generation
- start values

If two implementations are reading different data, the benchmark is weak before it starts.

### Step 3: Run the production smoke benchmark

Before deep parity work, prove that the production JWAS path actually runs and produces expected outputs.

This step should verify:

- output files exist
- result tables contain expected columns
- unsupported configurations fail clearly

### Step 4: Align deterministic setup

Before comparing chains, prove:

- same priors
- same fixed or estimated settings
- same initial state
- same data after preprocessing

If needed, export:

- `config.csv`
- `initial_state.csv`
- `initial_scalars.csv`

### Step 5: Run a one-step controlled replay

If a reference implementation exists, do not jump directly from setup to long-chain parity.

First, isolate one iteration:

- same initial state
- same draws
- same substeps written to files

This should answer whether the update math itself is aligned.

### Step 6: Compare production JWAS against local Julia diagnostics

If the one-step replay looks correct, compare:

- production `runMCMC(...)`
- local Julia driver for the same method

This tells you whether any remaining mismatch is internal to JWAS or only cross-language.

### Step 7: Run long-chain parity

Only after the earlier layers are clean should you use long chains to compare posterior summaries.

Long-chain parity should be the basis for the final correctness statement.

### Step 8: Run multiseed parity

If the method is stochastic and important, this should be the default final step.

Multiseed parity answers:

- is the observed difference stable?
- is it larger than ordinary between-chain variability?

### Step 9: Write the benchmark report

Every benchmark series should end with a short markdown note that states:

- what was benchmarked
- how it was benchmarked
- what matched
- what did not match
- whether the remaining differences are acceptable

## Required Artifacts

The minimum set of artifacts for a serious JWAS prototype-to-production benchmark is:

- benchmark input data or data-generation recipe
- saved config and starting state
- scripts used to run JWAS
- scripts used to run any reference implementation
- comparison script
- scalar summary tables
- method-specific tables such as `pi`
- markdown report

Recommended directory split:

- benchmark drivers and helpers: `benchmarks/`
- concrete benchmark conclusions: `benchmarks/reports/`
- manual guidance: `docs/src/manual/`

## Common Failure Modes

The following are common mistakes and should be treated as red flags.

### Comparing raw stochastic traces across languages without shared draws

This usually tells you more about different RNG streams than about algorithm correctness.

### Concluding “bug” from one short chain

For Bayesian methods, a single short chain is often not enough evidence.

### Benchmarking only a local helper instead of production `runMCMC`

This can miss integration problems in validation, initialization, posterior accumulation, or output code.

### Comparing only one metric

A single scalar such as `sigmaSq` is not enough. At minimum, compare:

- marker variance scale
- residual variance
- sparsity or nonzero frequency
- method-specific outputs such as `pi`

### Changing preprocessing between implementations

If marker filtering, centering, ordering, or missing-value handling differ, parity conclusions are weak.

### Mixing correctness and performance questions

First show the method is correct. Then measure speed and memory.

## Acceptance Checklist

Before saying a new method or backend is benchmarked for JWAS, confirm:

- the production JWAS path was benchmarked
- data and priors were aligned across implementations
- a smoke benchmark passed
- deterministic setup was checked
- a one-step replay was used when a reference implementation existed
- long-chain parity was run
- multiseed parity was considered for final conclusions
- outputs and benchmark scripts were saved
- a markdown benchmark note was written

If any of these are missing, the benchmark may still be useful for debugging, but it should not be treated as the final validation of production behavior.

## Worked Example: BayesR

BayesR provides a good example of why JWAS benchmarking must be layered.

The BayesR work started from:

- an external R reference script
- a new JWAS production implementation
- a parity question focused on posterior `sigmaSq`, residual variance, sparsity, and `pi`

The useful sequence was:

1. run production JWAS and the R reference on the same exported dataset
2. observe that a short single-chain parity run showed a visible `sigmaSq` mismatch
3. align priors and starting values explicitly
4. run a controlled one-step replay with shared draws
5. show that the one-step BayesR math matched closely
6. compare production JWAS against benchmark-local Julia and show those matched closely
7. test whether the short-chain mismatch was just Monte Carlo noise
8. switch to long-chain multiseed parity

Main lesson:

- the early mismatch looked like a bug
- the layered benchmark process showed it was mostly a short-chain Monte Carlo effect

That is exactly why prototype-to-production benchmarking in JWAS must not stop at one short stochastic comparison.

## Practical Recommendation

When in doubt, follow this rule:

- use short benchmarks to debug
- use controlled replay to validate formulas
- use production long-chain multiseed summaries to make final claims

This is the standard JWAS path for turning prototype code into production package functionality.
