# BayesR3 Block Strategy in JWAS (`fast_blocks`)

JWAS adopts and generalizes the BayesR3 computational strategy for block marker sampling. On this page, that strategy means three operations: constructing a block right-hand side (RHS), performing repeated conditional marker updates through the block Gram matrix, and reconciling the global corrected phenotype / residual when leaving the block.

This terminology describes the shared computational strategy; JWAS does not claim to reproduce every statistical-model or scheduling detail of the original BayesR3 software. JWAS uses the strategy in the supported single-trait BayesA/B/C and BayesR paths, annotated variants, and multi-trait paths listed below.

The default execution schedule remains an exact sequential block sweep. The explicit `independent_blocks=true` option instead selects an approximate parallel schedule.
For detailed non-block vs block memory accounting, see [Memory Usage](memory_usage.md).
For a real cluster timing benchmark at `N=50,000` targeting `P=2,000,000` and `chain_length=2000`, see [Benchmark](benchmark.md).

[BayesR3 paper (Methods)](https://www.nature.com/articles/s42003-022-03624-1)

## Four Independent Choices

- **Statistical model:** BayesA/B/C, BayesR, or an annotated or multi-trait variant determines the prior and conditional marker-sampling logic.
- **Block partition:** `fast_blocks` determines which markers are grouped for the block linear algebra. JWAS can choose a size, use a requested fixed size, or use explicit contiguous blocks.
- **Repeat policy:** the sampler determines how many conditional marker-update passes are performed while a block RHS and Gram matrix are active.
- **Execution schedule:** `independent_blocks=false` updates blocks sequentially and reconciles the global residual after each block; `independent_blocks=true` updates blocks from a frozen sweep-level residual and merges their changes after a barrier.

Changing the statistical model from BayesR to BayesC changes its prior and conditional sampling details, but not the core block RHS, Gram-matrix update, and block-exit residual-reconciliation strategy. Likewise, changing the block partition does not by itself select the independent-block approximation; only `independent_blocks=true` changes the inter-block execution schedule.

## User-Defined Unequal Blocks

Users may choose block boundaries using scientific information about LD, recombination, IBD, pedigree information, or genomic regions. These user-defined blocks do not need to have equal sizes.

The current API accepts ordered block start positions, so every block is contiguous in the current marker order. Arbitrary non-contiguous marker membership is not currently supported.

Grouping strongly coupled markers in the same block can improve computational locality or mixing per unit of computation. This does not guarantee better convergence: mixing and convergence depend on posterior correlations, block sizes, boundaries, and the repeat policy.

## Configuration and Supported Paths

Block updates are enabled with:

```julia
out = runMCMC(model, phenotypes; fast_blocks=true)
# or provide a fixed block size
out = runMCMC(model, phenotypes; fast_blocks=64)
# or provide explicit marker block starts
out = runMCMC(model, phenotypes; fast_blocks=[1, 501, 975, 1600])
# or opt into the approximate independent-block mode
out = runMCMC(model, phenotypes; fast_blocks=[1, 501, 975, 1600],
              independent_blocks=true)
```

### Quick Semantics Checklist

- `fast_blocks=true` chooses `block_size = floor(sqrt(nObs))`; numeric `fast_blocks=<block_size>` uses the requested fixed size. In the current implementation, a numeric block size should be less than `nMarkers` because chain-length scaling indexes the second block start.
- An explicit `fast_blocks=[...]` vector contains ordered marker start positions, not block lengths. The positions must be sorted, unique, start at `1`, and stay within `1:nMarkers`. In the current marker order, each position starts a contiguous block, so the vector can define unequal block sizes. Explicit starts use full-sweep semantics: one outer MCMC iteration sweeps every supplied block, and JWAS does not rescale `chain_length`.
- `independent_blocks=false` is the default exact sequential schedule: JWAS refreshes the global corrected phenotype / residual after each block. `independent_blocks=true` instead freezes the sweep-level residual, updates blocks with Julia threads, and merges their deltas after a barrier; this parallel schedule is approximate unless off-block weighted genotype crossproducts (`X_b' W X_c` for `b != c`) are effectively zero.
- `fast_blocks` is currently dense-storage only; `storage=:stream` rejects `fast_blocks != false`.
- On servers, set `OPENBLAS_NUM_THREADS=1` when using Julia threads (for example, with `independent_blocks=true`) to avoid nested BLAS oversubscription.

### Supported Sampler Paths

- In single-trait BayesA/B/C, JWAS calls `BayesABC_block!`.
- In single-trait BayesR, JWAS calls `BayesR_block!`.
- Dense annotated single-trait BayesC and BayesR use the same block machinery with marker-specific annotation priors.
- In multi-trait BayesA/B/C with unconstrained marker covariance (`Mi.G.constraint == false`), JWAS calls `MTBayesABC_block!` and honors `multi_trait_sampler=:I`, `:II`, or `:auto`.
- Dense annotated 2-trait BayesC uses the same multi-trait block sampler choices.
- If `Mi.G.constraint == true`, JWAS uses `megaBayesABC!` (non-block path).

## Single-Trait BayesC Example: Without Blocks

In the standard BayesC update (`BayesABC!`), each marker `j` is updated one-by-one:

1. Compute conditional posterior terms (`rhs`, `lhs`, `gHat`) from current `yCorr`.
2. Compute marker inclusion probability `probDelta1`.
3. Sample `δ[j]` (include/exclude) and sample/update `β[j]`, `α[j]`.
4. Immediately update `yCorr` using marker column `x_j`.

So `yCorr` is updated every marker.

## Single-Trait BayesC Example: Block Strategy in JWAS

In the block version (`BayesABC_block!`), markers are partitioned into blocks.

### Precomputation (once before MCMC)

For each block `b`, JWAS builds:

- `X_b` (block genotype matrix)
- `X_b'R^{-1}X_b` (`XpRinvX`, block Gram matrix)

`X_b'R^{-1}` is not persisted as a block matrix. Instead, JWAS computes
`X_b'R^{-1}yCorr` on demand into a reusable block workspace vector each outer
iteration.

### Update flow inside each MCMC outer iteration

For each block:

1. Build block RHS once: `XpRinvycorr = X_b'R^{-1} yCorr`.
2. Update markers inside the block using BayesC logic.
3. Instead of touching full `yCorr` each marker, update the block RHS using columns of `X_b'R^{-1}X_b`.
4. After finishing the block, update full `yCorr` once:
   `yCorr += X_b * (α_old_block - α_new_block)`.

### Important implementation detail

In the current JWAS code:

- for `fast_blocks=true` or numeric `fast_blocks`, inner repeats are set to `nreps = block_size`;
- for `fast_blocks=true` or numeric `fast_blocks`, outer `chain_length` is reduced by approximately `block_size`;
- for explicit block starts, one MCMC iteration means one full sweep over the supplied blocks and JWAS does not rescale `chain_length`.

This keeps effective marker-update work on a similar scale while moving much of the per-marker work from `nObs`-length operations to `block_size`-length operations.
The explicit block-start form uses sweep semantics because user-provided LD or pedigree-informed blocks can have different sizes.

## Exact Sequential Blocks vs Independent Blocks

The default block sampler composes conditional Gibbs updates in sequence. After
one block finishes, JWAS reconciles that block's effect changes into the global
corrected phenotype / residual before it constructs the right-hand side for the
next block. These steps form the exact sequential blocked schedule for the
implemented model.

This is a different transition schedule from dense marker-by-marker Gibbs: it
groups conditional marker updates by block and can repeat updates within a
block. It does not, however, introduce the frozen-residual approximation
described below. That approximation is selected only by
`independent_blocks=true`.

### Single-Trait View

Let:

- `y*` be the phenotype after non-marker terms are removed
- `X = [X_1, X_2, ..., X_B]` be the marker matrix partitioned into blocks
- `alpha_b` be marker effects in block `b`
- `W` be the diagonal observation-weight matrix represented by `Rinv`

The sequential update for block `b` uses:

```text
s_b = X_b' W (r + X_b alpha_b)
    = X_b' W (y* - sum_{c != b} X_c alpha_c)
```

By contrast, `independent_blocks=true` freezes the global residual at the start
of a marker sweep, so blocks do not see effect changes made by other blocks
during that sweep. This independent-block schedule is exact only when the
off-block weighted genotype crossproducts vanish:

```text
X_b' W X_c = 0  for b != c
```

Under that condition, the off-block terms vanish. For two marker blocks, this is the user's intuition:

```text
x1' W (y - x1 alpha1 - x2 alpha2) = x1' W (y - x1 alpha1)
```

when:

```text
x1' W x2 = 0
```

If the off-block crossproducts are small but not zero,
`independent_blocks=true` is an approximation.

### Multi-Trait View

For multi-trait models, let:

- `Y*` be the matrix of corrected phenotypes after non-marker terms are removed
- `A_b` be the marker-effect matrix for block `b`
- `R^{-1}` be the residual precision matrix among traits

The within-block sampler still uses the existing multi-trait residual covariance logic through `R^{-1}`.
The independent-block assumption is still about genotype-side block leakage:

```text
X_b' W X_c = 0  for b != c
```

So the independent-block approximation:

1. freezes the trait-wise corrected phenotype vectors at the start of the marker sweep;
2. updates each genotype block independently from that frozen state;
3. applies all block effect deltas to the global corrected phenotypes after all blocks finish.

The prior is not changed by this option. Plain and annotated priors keep their existing marker-state logic.

### Server Use

The independent-block path uses Julia threads over blocks.
On a server, use:

```bash
export JULIA_NUM_THREADS=<num_cores>
export OPENBLAS_NUM_THREADS=1
julia --project=.
```

`OPENBLAS_NUM_THREADS=1` avoids oversubscribing the machine with nested BLAS threads while JWAS is already parallelizing over marker blocks.

### Detailed Comparison with BayesR3

Reference paper: BayesR3 (Communications Biology, 2022), DOI: `10.1038/s42003-022-03624-1`.

The tables below are a worked comparison between the original BayesR3 method and the JWAS block BayesC implementation. JWAS also applies this generalized block strategy to BayesR and to the supported annotated paths described above.

For the repeat schedules in this comparison:

- `b` is the nominal block size.
- `s_i` is the realized number of markers in block `i`.
- `r_i` is the number of inner marker-update repeats performed while block `i` is active.

#### Step-by-step correspondence

| BayesR3 paper step | JWAS block BayesC implementation | Status |
| --- | --- | --- |
| Partition markers into blocks | `fast_blocks` builds marker blocks | Same strategy |
| Build per-block RHS \(r_b = V_b'We\) | `block_rhs!(XpRinvycorr, XArray[i], yCorr, Rinv, unit_weights)` | Same strategy |
| Within-block marker update uses current block RHS | BayesC per-marker `rhs/lhs/gHat` from `XpRinvycorr` | Same strategy |
| In-block RHS correction via block Gram column | `BLAS.axpy!(..., view(XpRinvX[i],:,j), XpRinvycorr)` | Same strategy |
| Update residual once on block exit | `yCorr += X_b*(α_old_block-α_new_block)` | Same strategy |

#### What is different

| Topic | BayesR3 paper | JWAS block BayesC | Practical implication |
| --- | --- | --- | --- |
| Marker prior model | BayesR mixture (multiple non-zero normal components plus zero component) | BayesC spike-slab style inclusion (`δ∈{0,1}` for this path) | Same acceleration idea, different posterior model |
| Marker state sampling | Multi-class mixture state | Binary include/exclude state | Not numerically identical to BayesR |
| Inner-repeat schedule | `r_i = b` for every block, including a final block with `s_i < b` | `r_i = s_i` | Short blocks receive fewer repeats; explicit unequal blocks can receive different repeat counts |
| Outer-loop scheduling | Described as a fixed block sweep schedule | For `fast_blocks=true` or numeric `fast_blocks=<block_size>`, JWAS rescales outer `chain_length` by block size; explicit starts keep full-sweep `chain_length` | Compare effective updates and partition semantics, not just outer iterations |
| Scope | BayesR algorithm | JWAS block path is wired to BayesA/B/C and BayesR marker samplers, with multi-trait BayesC honoring sampler I/II dispatch | Strategy reused across Bayesian alphabet members |

#### Scheduling detail (explicit)

Under the original BayesR3 convention, most blocks have `s_i = b`, while the final block may have `s_i < b`. BayesR3 nevertheless uses `r_i = b` for that final short block. Consequently, every SNP effect receives the same number of updates under the paper's chain-length convention.

Current JWAS repeat policies depend on the marker sampler:

- BayesC uses `r_i = s_i`.
- BayesR uses `r_i = 1` during burn-in and `r_i = s_i` afterward.

With fixed-size partitioning, full blocks have `s_i = b`, while a final short block has a smaller repeat count under the JWAS BayesC policy and under the post-burn-in JWAS BayesR policy. Explicit unequal starts can likewise produce different repeat counts among blocks. These scheduling differences do not stop JWAS from using the generalized BayesR3 block strategy: the block RHS, Gram-matrix update, and block-exit residual reconciliation remain the shared computational pattern.

## Algorithm Comparison

| Aspect | Standard BayesC (`BayesABC!`) | Exact block BayesC (`BayesABC_block!`) | Independent blocks (`independent_blocks=true`) |
| --- | --- | --- | --- |
| Update unit | One marker at a time | One block, then markers inside block | One block per thread/chunk, then markers inside block |
| `yCorr` / corrected phenotype updates | Every marker | Once per block, before the next block | Once after all blocks in the sweep finish |
| Inter-block schedule | Fully sequential | Sequential and exact | Parallelizable and approximate unless off-block crossproducts vanish |
| Main per-marker linear algebra size | `nObs` | `block_size` (inside block cache) | `block_size` (inside block cache) |
| Extra precompute | Minimal | `X_b'R^{-1}X_b` for all blocks (RHS computed on demand) | Same block precompute |
| Extra memory | Minimal | Stores block Gram matrices (`XpRinvX`) and block workspaces | Same plus block-local deltas and thread-local buffers |
| Chain behavior | Direct `chain_length` | Fixed-size blocks use inner repeats + outer chain scaling; explicit starts use full sweeps | Same chain semantics as the selected `fast_blocks` partition |

## Computational Complexity

Use the notation:

- `N`: number of records (`nObs`)
- `P`: number of markers (`nMarkers`)
- `b`: nominal block size
- `B = ceil(P/b)`: number of blocks
- `L`: requested chain length

### Standard BayesC (non-block)

- Per MCMC iteration: `O(PN)` (marker-wise dot products and residual updates over `N` records)
- Total over `L` iterations: `O(LPN)`

### JWAS block BayesC (current implementation)

Let block sizes be `s_i` with `sum_i s_i = P`.

Per outer iteration:

- Block RHS construction across all blocks: `O(NP)`
- In-block updates: `O(sum_i s_i^3)` (because `nreps = s_i` and in-block RHS updates are length-`s_i`)
- Residual updates on block exit across all blocks: `O(NP)`

So per outer iteration:

- `O(NP + sum_i s_i^3)`
- With near-uniform blocks (`s_i ≈ b`): `O(NP + P b^2)`

For `fast_blocks=true` or numeric `fast_blocks=<block_size>`, JWAS rescales outer iterations to approximately `m = floor(L/b)`, so the main total cost is:

- `O((L/b) * (NP + P b^2)) = O(LP(N/b + b))`

For explicit `fast_blocks=[...]`, JWAS instead treats one outer iteration as one full sweep over the supplied block starts and does not rescale `chain_length`. The total over `L` requested full sweeps is:

- `O(L * (NP + sum_i s_i^3))`

### JWAS arithmetic operation-count example (`N=200,000`, `P=2,000,000`)

Assume `fast_blocks=true`, so JWAS uses `b = floor(sqrt(N)) = 447`.
Then:

- `B = ceil(P/b) = ceil(2,000,000/447) = 4,475` blocks
- Standard BayesC total scaling: `O(LPN) = O(L * 2,000,000 * 200,000)`
- JWAS block BayesC main scaling: `O(LP(N/b + b)) = O(L * 2,000,000 * (200,000/447 + 447))`

Dividing these leading arithmetic expressions gives
`N / (N/b + b) ≈ 224`. This is a theoretical arithmetic operation-count
comparison between the two JWAS paths, not a measured speedup. Wall-clock
performance also depends on implementation constants, memory traffic, BLAS,
hardware, and the workload.

### BayesR3 published empirical runtime fit

BayesR3 uses the same block-strategy family but is a distinct BayesR mixture
implementation with its own repeat schedule. Figure 5 of the BayesR3 paper
reports that its measured processing time per SNP is proportional to
`(N + b)/b = N/b + 1`. This is an empirical fit for that implementation, data,
software libraries and hardware, and the sample-size and block-size regime
tested in the paper. It is not an operation-count derivation.

The JWAS operation-count expression above and the BayesR3 empirical fit must
not be compared coefficient by coefficient or interpreted as a JWAS-versus-
BayesR3 software benchmark. A software performance comparison would require
both implementations to be timed under a controlled benchmark using aligned
models, chain schedules, data, libraries, and hardware.

## Detailed Resource Model (Current `fast_blocks` Path)

Use:

- `N`: records
- `P`: markers
- `b`: nominal block size
- `s_i`: size of block `i`, `sum_i s_i = P`
- `t`: bytes per floating-point value (`4` for `Float32`, `8` for `Float64`)

### Memory formulas

Current block implementation (after removing persistent `XRinvArray`) stores:

- dense `X`: `N*P` values
- block Gram matrices `XpRinvX`: `sum_i s_i^2` values
- marker summary `xpRinvx`: `P` values
- optional `xRinvArray` extra `N*P` for non-unit weights
- small reusable block workspaces (`O(b)`) for RHS and local deltas

Approximate totals:

- Unit weights:
  - `Mem_block_unit ~= t * (N*P + sum_i s_i^2 + P) + O(b*t)`
- Non-unit weights:
  - `Mem_block_nonunit ~= t * (2*N*P + sum_i s_i^2 + P) + O(b*t)`

### Runtime working set

Per block update, active high-volume buffers are:

- `X_b` view (`N x s_i`; no data copy)
- `XpRinvX[i]` (`s_i x s_i`)
- block RHS workspace (`s_i`)
- `yCorr` (`N`)

So peak additional block-local workspace is roughly:

- `O(s_i^2 + s_i + N)` values

### I/O and precompute considerations

For in-memory dense `X`, no out-of-core read is required during MCMC sweeps.
The expensive setup term is building `XpRinvX`:

- precompute cost scales with approximately `O(N * sum_i s_i^2)` (`~O(NPb)` under near-uniform blocks)

This setup can dominate startup time for very large `N,P`, even when per-iteration sampling is fast.

## Worked Large-Scale Example (`N=500,000`, `P=2,000,000`)

Assume `fast_blocks=true`, so `b=floor(sqrt(N))=707`.

- `B = ceil(P/b) = 2,829`
- `sum_i s_i^2 = 1,413,937,788`

Memory-relevant terms:

| Term | Float32 | Float64 |
| --- | ---: | ---: |
| `X` | `~4.00 TB` (`3.64 TiB`) | `~8.00 TB` (`7.28 TiB`) |
| `XpRinvX` | `~5.66 GB` (`5.27 GiB`) | `~11.31 GB` (`10.53 GiB`) |
| `xpRinvx` | `~8.0 MB` | `~16.0 MB` |

So for unit weights, block mode remains dominated by `X`, with `XpRinvX` as the main incremental memory term.

## What To Watch Closely

1. `block_size` choice:
   - too small: less speedup (`N/b` term remains large)
   - too large: `XpRinvX` memory and precompute cost rise (`~P*b` memory and `~NPb` setup trend)
2. Effective chain length:
   - with `fast_blocks=true` or numeric `fast_blocks=<block_size>`, current implementation rescales outer iterations and uses inner repeats (`nreps = block_size`)
   - with explicit `fast_blocks=[...]`, one outer iteration is a full sweep over the supplied block starts and `chain_length` is not rescaled
   - compare runs by effective updates, not only outer-iteration count
3. Final short block behavior:
   - last block uses smaller `nreps` (equal to its own size), so sweep symmetry differs slightly
4. Independent blocks:
   - exact only when off-block weighted genotype crossproducts vanish
   - otherwise it is an explicit approximation for speed and parallelism
5. Multi-trait block path specifics:
   - sampler I, sampler II, and `:auto` are supported for unconstrained covariance mode
   - extra temporaries (e.g., trait-wise old-alpha handling) can become noticeable at larger trait counts
6. Numerical reproducibility:
   - mathematically equivalent refactors (e.g., in-place BLAS updates) can change floating-point roundoff
   - expect tiny non-bitwise differences, especially in `Float32`
7. Weighting mode:
   - block `XRinvArray` is no longer persisted
   - separate non-unit weighted `xRinvArray` materialization remains a distinct issue
8. Scope constraints:
   - current source notes this fast block option is intended for one genotype category
   - numeric `fast_blocks` should keep `block_size < nMarkers`

## Example: Speed/Memory Tradeoff

Assume:

- `nObs = 5_000`
- `nMarkers = 50_000`
- `block_size = 70` (similar to `sqrt(nObs)`)
- `nBlocks = 715`

Approximate memory (Float32):

| Item | Approx size |
| --- | --- |
| Genotype matrix `X` (`nObs x nMarkers`) | `~953.7 MiB` |
| Extra `XpRinvX` in block mode | `~13.4 MiB` |

Interpretation:

- Block mode can be much faster for large `nObs`, because heavy per-marker operations are shifted to block-level cached operations.
- Block mode uses more memory than non-block mode, mainly from `XpRinvX` (plus small block workspaces).
- Practical speedup is typically below theoretical arithmetic speedup due to random branching, allocation overhead, and BLAS/runtime effects.

## Practical Guidance

1. Start with `fast_blocks=true` for large marker sets and enough RAM.
2. If memory is tight, set a smaller numeric block size (e.g., `32` or `64`) and benchmark.
3. If speed gain is small, try a few block sizes and choose based on wall time + memory headroom.
4. If you have external LD, recombination, IBD, or pedigree-informed blocks, pass them as explicit starts with `fast_blocks=[...]`.
5. Use `independent_blocks=true` only when you intentionally accept the independent-block approximation and want block-level thread parallelism.
