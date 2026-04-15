# Jacobi Fast-Block Parallelism Design

**Status:** proposed

**Goal:** add an explicit approximate fast-block mode that assumes user-supplied genotype blocks are independent enough to support Jacobi-style block updates and server-side parallel execution.

## Motivation

Current JWAS fast-block sampling is a **Gauss-Seidel-style** block sweep:

- sample block `1`
- update the global corrected phenotype / residual state
- sample block `2` using that updated state
- continue sequentially

This is exact for the current model, but it is inherently serial across blocks. The main performance bottleneck is therefore not the within-block algebra; it is the fact that blocks cannot be updated simultaneously.

The new feature is an explicit **Jacobi-style** alternative:

- freeze the sweep-level corrected phenotype / residual state at the beginning of the sweep
- sample every block against that same frozen snapshot
- merge the block updates at the end of the sweep

This makes block-parallel execution natural on a multi-core server. It is exact only when the supplied blocks are conditionally independent in the relevant crossproduct structure; otherwise it is an explicit approximation.

## Scope

This design covers all current fast-block BayesianAlphabet paths:

- single-trait BayesC / BayesB / BayesA via `BayesABC_block!`
- single-trait annotated BayesC via the same path
- single-trait BayesR via `BayesR_block!`
- single-trait annotated BayesR via the same path
- multi-trait BayesC via `MTBayesABC_block!`
- annotated 2-trait BayesC via the same path

This design does **not** add:

- automatic block construction
- pedigree-driven block construction inside JWAS
- distributed multi-process execution
- multi-trait BayesR support

For v1, JWAS only consumes a user-provided block partition or a fixed-width partition derived from existing `fast_blocks` behavior.

## Statistical Framing

### Single-trait

Let:

- `y*` be the phenotype after all non-marker terms are removed
- `X = [X_1, X_2, ..., X_B]` be the marker matrix partitioned into blocks
- `alpha = (alpha_1', ..., alpha_B')'`
- `W` be the diagonal observation-weight matrix currently encoded by `Rinv` in the single-trait block code

The full single-trait residual is:

- `r = y* - sum_b X_b alpha_b`

The block-local weighted score for block `b` is:

- `s_b = X_b' W (r + X_b alpha_b)`

Equivalently,

- `s_b = X_b' W (y* - sum_{c != b} X_c alpha_c)`

#### Gauss-Seidel block sweep

Current JWAS fast blocks are Gauss-Seidel-style:

- when updating block `b`, previously updated blocks in the same sweep already change the global residual
- the current block therefore sees the newest available values of earlier blocks

This is the exact current sampler.

#### Jacobi block sweep

The proposed mode freezes `alpha_old` and the corresponding sweep residual:

- `r_old = y* - sum_b X_b alpha_b_old`

Then every block uses:

- `s_b^J = X_b' W (r_old + X_b alpha_b_old)`

All blocks are updated against the same frozen sweep snapshot.

This becomes exact if the off-block weighted crossproducts vanish:

- `X_b' W X_c = 0` for all `b != c`

Otherwise the method is an explicit approximation.

### Multi-trait

Let:

- `Y*` denote the corrected trait matrix after non-marker terms are removed
- `X = [X_1, ..., X_B]` be the genotype matrix partitioned into blocks
- `A_b` be the marker-effect matrix for block `b`
- `R` be the residual covariance across traits

The current multi-trait dense and block samplers use the trait coupling through `R^{-1}` and the genotype-side weighted crossproducts `x' R^{-1} x` / `X' R^{-1} X`.

The relevant off-block independence statement is therefore:

- the off-block contribution to the full precision term `X' R^{-1} X` vanishes

Operationally, for a given genotype partition, this means the off-block genotype crossproducts should be negligible after weighting:

- `X_b' W X_c ≈ 0` for `b != c`

where `W` is the observation-weight structure induced by the current implementation.

#### Gauss-Seidel block sweep

Current `MTBayesABC_block!` is Gauss-Seidel-style:

- block-local RHS objects are built from the current global corrected phenotypes
- after each block finishes, the global trait-wise corrected phenotype vectors are updated
- later blocks therefore depend on earlier blocks within the same sweep

This is the exact current sampler.

#### Jacobi block sweep

The proposed multi-trait Jacobi mode freezes the sweep-level corrected phenotype vectors at the start of the sweep and updates all blocks against that same snapshot.

Within each block:

- the current multi-trait BayesC machinery remains unchanged
- sampler I vs sampler II remains unchanged
- plain vs annotated priors remain unchanged

What changes is only the inter-block coupling:

- no within-sweep global residual correction between blocks
- one reconciliation step after all blocks finish

If the off-block contributions to `X' R^{-1} X` vanish, this is exact under the independence assumption. Otherwise it is an approximate Jacobi-style sampler.

## Computational Framing

### Current fast-block mode

Current fast blocks are exact but serial across blocks:

1. build block-local RHS from the current corrected phenotype / residual
2. sample one block
3. update the global corrected phenotype / residual
4. continue to the next block

This is a textbook Gauss-Seidel update order.

### Proposed Jacobi mode

Jacobi mode changes the sweep order:

1. compute the frozen sweep snapshot
2. launch independent block updates from that same snapshot
3. collect each block's effect delta
4. apply the combined delta to the global corrected phenotype / residual at the barrier

This is the natural form for server-side block parallelism.

### Server parallelism

For a single multi-core server, the intended execution model is Julia threads:

- one Julia process
- many Julia threads
- BLAS threads pinned to `1`
- one barrier per MCMC sweep

Recommended server environment:

```bash
export JULIA_NUM_THREADS=<num_cores>
export OPENBLAS_NUM_THREADS=1
julia --project=.
```

Each worker thread should own:

- one block or a chunk of blocks
- thread-local RHS buffers
- thread-local temporary arrays
- thread-local RNG

Shared writes should be avoided during the block update itself. Each block should return its local effect delta, and the barrier reduction should apply those deltas once per sweep.

### Expected performance relative to OpenMP

This Jacobi formulation is algorithmically very similar to what one would want in an OpenMP implementation.

Expected behavior:

- the major speedup comes from changing **serial Gauss-Seidel blocks** into **parallel Jacobi blocks**
- Julia threads can be in the same performance ballpark as OpenMP if:
  - blocks are coarse enough
  - thread-local allocations are controlled
  - BLAS oversubscription is avoided
  - synchronization is limited to one barrier per sweep

What will still matter:

- memory bandwidth
- block-size imbalance
- temporary allocations
- thread scheduling overhead

So the honest performance statement is:

- same order of magnitude as OpenMP is plausible
- parity is not guaranteed
- the algorithmic change matters more than the language choice

## API Design

### User-facing options

Add an explicit block sweep option to `runMCMC(...)`:

- `block_sweep::Symbol = :gauss_seidel`

Accepted values:

- `:gauss_seidel` = exact current fast-block sweep
- `:jacobi` = explicit approximate independent-block sweep

This option only matters when `fast_blocks != false`.

### Block partition input

Extend `fast_blocks` so the user can provide the partition directly.

Accepted forms:

- `false` = no block mode
- `true` = existing heuristic block size based on `sqrt(nObs)`
- `Number` = existing fixed block size
- `AbstractVector{<:Integer}` = explicit block-start vector supplied by the user

The vector form should be interpreted exactly as current internal `block_starts`:

- sorted
- unique
- 1-based
- first entry must be `1`
- entries must be within `[1, nMarkers]`

This keeps the new feature aligned with the user's stated workflow: JWAS consumes a block partition provided by the user.

## Pedigree-Informed Blocking

JWAS will **not** build pedigree-derived blocks in v1.

However, pedigree can still help users construct better supplied blocks outside JWAS:

- small effective population size implies long shared haplotypes
- pedigree and family structure can help avoid splitting inside long co-segregating regions
- identity-by-descent or family-aware recombination information can suggest better contiguous breakpoints

The important scientific distinction is:

- pedigree may help propose boundaries
- the sampler's actual independence criterion is still block leakage in genotype crossproducts

So JWAS should document pedigree-informed partitioning as user guidance, not as built-in functionality.

## Quality of the Independence Assumption

For this Jacobi mode, “better blocks” means:

- large within-block dependence
- small off-block dependence

The relevant diagnostics are based on off-block weighted crossproducts, not just casual LD language.

Useful leakage summaries include:

- max off-block absolute crossproduct
- sum of absolute off-block crossproducts
- adjacent-boundary leakage
- ratio of off-block mass to within-block mass

These diagnostics should be exposed in documentation and, if feasible, via helper output in a later iteration. They are not required for the first implementation.

## Method Coverage

The same Jacobi/Gauss-Seidel distinction should be implemented consistently for:

- single-trait BayesABC block samplers
- single-trait BayesR block samplers
- multi-trait MTBayesABC block samplers

The prior logic should remain unchanged:

- plain vs annotated BayesC uses the same prior abstraction already present
- BayesR vs annotated BayesR keeps its existing prior-row logic
- sampler I vs sampler II logic for multi-trait BayesC remains unchanged inside each block

Only the inter-block sweep order changes.

## Verification Strategy

The design needs two kinds of evidence.

### Exactness under block independence

On synthetic block-diagonal test data where off-block crossproducts are exactly zero:

- Gauss-Seidel and Jacobi should agree up to Monte Carlo error
- single-trait and multi-trait paths should both satisfy this

### Approximation behavior under realistic coupling

On realistic livestock-style data with nonzero cross-block coupling:

- Jacobi should be benchmarked as an approximate method
- results should be compared to Gauss-Seidel on:
  - posterior summaries
  - EBV / prediction accuracy
  - runtime
  - strong scaling across server threads

## Recommended v1

1. Add `block_sweep=:gauss_seidel|:jacobi` to `runMCMC`.
2. Allow explicit user-provided block-start vectors in `fast_blocks`.
3. Implement Jacobi sweeps for:
   - single-trait BayesABC
   - single-trait BayesR
   - multi-trait MTBayesABC
4. Use Julia threads for within-sweep block parallelism.
5. Treat Jacobi as an explicit approximate option in the manual and benchmarks.

This is the smallest design that is technically honest and directly useful for fast server parallel computing.
