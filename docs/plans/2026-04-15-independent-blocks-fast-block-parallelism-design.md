# Independent-Blocks Fast-Block Parallelism Design

**Status:** proposed

**Goal:** add an explicit approximate fast-block mode, `independent_blocks=true`, that assumes user-supplied marker blocks are independent enough to update in parallel within each MCMC sweep.

## Motivation

Current JWAS fast-block sampling is an exact sequential block sweep:

- sample block `1`
- update the global corrected phenotype / residual state
- sample block `2` using that updated state
- continue one block at a time

This is statistically exact for the current model, but it is serial across blocks.

The proposed `independent_blocks=true` mode changes the inter-block update:

- freeze the corrected phenotype / residual state at the start of the marker sweep
- update each block from that same frozen snapshot
- merge all block effect changes at the end of the sweep

This makes the blocks parallelizable on a multi-core server.

The mode is exact only when the supplied blocks are independent in the relevant weighted genotype crossproduct. Otherwise it is an explicit approximation.

## User-Facing API

Add a new `runMCMC` keyword:

```julia
independent_blocks::Bool = false
```

Behavior:

- `independent_blocks=false`: current exact fast-block behavior
- `independent_blocks=true`: approximate independent-block sweep
- `independent_blocks=true` requires `fast_blocks != false`

This name is intentionally practical. Users do not need to know numerical-analysis terminology. The option says exactly what assumption they are choosing.

## Block Partition Input

Current `fast_blocks` supports:

- `false`: no block mode
- `true`: JWAS chooses a heuristic block size
- `Number`: fixed block size

Extend `fast_blocks` to also accept:

- `AbstractVector{<:Integer}`: explicit block-start positions supplied by the user

The vector form should be validated as:

- sorted ascending
- unique
- first entry is `1`
- all entries are in `[1, nMarkers]`

Example:

```julia
runMCMC(model, data;
        fast_blocks=[1, 501, 975, 1600],
        independent_blocks=true)
```

This matches the intended livestock workflow: users can provide blocks chosen from LD boundaries, recombination patterns, pedigree-informed haplotype structure, or external software.

## Chain-Length Semantics

This needs to be explicit because current fixed-size fast blocks rewrite `chain_length`.

### Existing fixed-size fast blocks

For `fast_blocks=true` or `fast_blocks=<Number>`, preserve existing behavior:

- JWAS derives fixed-width block starts
- JWAS rescales the internal outer loop count as it does today
- existing users get the same behavior

### User-provided block starts

For `fast_blocks::AbstractVector{<:Integer}`, use sweep semantics:

- one MCMC iteration means one full sweep over all supplied blocks
- do **not** divide `chain_length` by a block size
- use one within-block pass per sweep unless a later feature adds an explicit block-repeat option

Reason:

- user-supplied LD blocks can have different sizes
- there is no single block size that makes the current `chain_length / block_size` rule correct
- sweep semantics are easier to explain and safer for variable-size blocks

This behavior should be documented clearly because it differs from the legacy fixed-size fast-block scaling.

## Statistical Framing

### Single-Trait

Let:

- `y*` be the phenotype after all non-marker terms are removed
- `X = [X_1, X_2, ..., X_B]` be the marker matrix partitioned into blocks
- `alpha = (alpha_1', ..., alpha_B')'`
- `W` be the diagonal observation-weight matrix represented by `Rinv` in the current block code

The full residual is:

```text
r = y* - sum_b X_b alpha_b
```

For block `b`, the current exact block update uses:

```text
s_b = X_b' W (r + X_b alpha_b)
```

This equals:

```text
s_b = X_b' W (y* - sum_{c != b} X_c alpha_c)
```

Under independent blocks:

```text
X_b' W X_c = 0  for all b != c
```

Then the off-block terms vanish. In the simple two-marker statement:

```text
x1' W (y - x1 alpha1 - x2 alpha2)
  = x1' W (y - x1 alpha1)
```

when:

```text
x1' W x2 = 0
```

That is the core assumption.

If the off-block weighted crossproducts are only small, not zero, then `independent_blocks=true` is an approximation.

### Multi-Trait

Let:

- `Y*` be the matrix of corrected phenotypes after non-marker terms are removed
- `X = [X_1, ..., X_B]` be the genotype matrix partitioned into blocks
- `A_b` be the marker-effect matrix for block `b`
- `R` be the residual covariance among traits

The within-block multi-trait conditional still uses the existing trait coupling through `R^{-1}`.

The inter-block independence condition is still driven by genotype-side block leakage:

```text
X_b' W X_c = 0  for all b != c
```

The shorthand “off-block pieces of `X'R^{-1}X` vanish” is useful conceptually, but the implementation should describe the actual quantity as weighted genotype block crossproducts. Trait coupling remains inside the block update through `R^{-1}`.

So `independent_blocks=true` means:

- freeze the trait-wise corrected phenotype vectors at the start of the marker sweep
- update each genotype block independently using the frozen vectors
- apply all block effect changes to the global corrected phenotypes at the end

The plain and annotated priors do not change this logic. The approximation is in the likelihood coupling across blocks, not in the prior.

## Parallel Computing Model

For a single server, the target implementation is shared-memory Julia threading.

Recommended runtime setup:

```bash
export JULIA_NUM_THREADS=<num_cores>
export OPENBLAS_NUM_THREADS=1
julia --project=.
```

Reason:

- parallelism should be over marker blocks
- BLAS should not also spawn many threads and oversubscribe the machine

Each thread should own:

- one block or a chunk of blocks
- local RHS buffers
- local temporary arrays
- local random number generator
- local block effect deltas

The threaded independent-block sweep should avoid shared writes during block updates. It should merge block deltas only at the end of the sweep.

## Performance Expectations

The main expected speedup is algorithmic:

- current exact fast blocks are sequential across blocks
- `independent_blocks=true` lets blocks run simultaneously

Compared with an OpenMP implementation:

- the same order of magnitude is plausible for coarse block work
- exact parity is not guaranteed
- block size, memory bandwidth, allocation behavior, and load balance matter more than the threading system alone

The practical path is:

1. implement the feature with Julia threads
2. benchmark strong scaling on the server
3. only consider lower-level OpenMP/C/Fortran work if profiling shows Julia threading overhead is the bottleneck

## Pedigree-Informed Blocking

JWAS should not construct pedigree-derived blocks in v1.

However, pedigree information can help users choose better external blocks:

- livestock populations often have small effective population size
- long haplotypes and family structure are common
- pedigree and identity-by-descent information can suggest boundaries that avoid splitting co-segregating regions

The scientific distinction is:

- pedigree can propose block boundaries
- the sampler’s independence assumption is still judged by genotype block leakage

So v1 should document this as guidance:

- use pedigree, IBD, recombination, LD, or external tools to propose blocks
- pass those blocks to JWAS through `fast_blocks=[...]`
- use benchmark diagnostics to evaluate whether the approximation is acceptable

## Block-Quality Diagnostics

Better independent blocks have:

- strong within-block dependence
- weak off-block dependence

Useful diagnostics include:

- max adjacent off-block weighted crossproduct
- sum of absolute off-block weighted crossproducts
- ratio of off-block mass to within-block mass
- comparison of independent-block posterior summaries to exact sequential fast-block summaries

These diagnostics are useful, but they are not required for v1 implementation. They can be added in a later benchmark or helper function.

## Method Coverage

The implementation should cover:

- single-trait BayesC / BayesB / BayesA
- single-trait annotated BayesC
- single-trait BayesR
- single-trait annotated BayesR
- multi-trait BayesC
- annotated 2-trait BayesC

This should be one conceptual feature:

- exact block mode: `independent_blocks=false`
- approximate parallel block mode: `independent_blocks=true`

The within-block sampler logic should remain unchanged.

## Verification Strategy

Use two classes of evidence.

### Block-Diagonal Synthetic Data

Construct toy data with exact block independence:

```text
X_b' W X_c = 0  for all b != c
```

On this data:

- exact sequential fast blocks and independent blocks should agree up to Monte Carlo error
- single-trait and multi-trait tests should both be covered

### Realistic Data

On realistic data, `independent_blocks=true` should be treated as approximate.

Benchmarks should compare:

- runtime
- server-thread scaling
- EBV / prediction accuracy
- marker-effect summaries
- posterior inclusion summaries

Marker-level comparisons must join by explicit marker keys, not by row order.

## Recommended v1

1. Add `independent_blocks::Bool=false` to `runMCMC`.
2. Allow explicit user-provided block starts in `fast_blocks`.
3. Preserve current exact fast-block behavior when `independent_blocks=false`.
4. Add independent-block sweeps for single-trait BayesABC, single-trait BayesR, and multi-trait MTBayesABC.
5. Implement independent blocks sequentially first.
6. Add thread-parallel execution after the sequential path is verified.
7. Benchmark approximation quality and server scaling.
