# Block BayesC (`fast_blocks`)

This page explains how JWAS implements block updates for BayesC marker sampling and how the block path changes speed and memory usage.
The block BayesC implementation uses a strategy similar to the blocked update scheme described in the BayesR3 paper.
For detailed non-block vs block memory accounting, see [Memory Usage](memory_usage.md).

BayesR3 paper (Methods): https://www.nature.com/articles/s42003-022-03624-1

## When This Path Is Used

Block updates are enabled with:

```julia
out = runMCMC(model, phenotypes; fast_blocks=true)
# or provide a fixed block size
out = runMCMC(model, phenotypes; fast_blocks=64)
```

- If `fast_blocks=true`, JWAS chooses `block_size = floor(sqrt(nObs))`.
- In single-trait BayesA/B/C, JWAS calls `BayesABC_block!`.
- In multi-trait BayesA/B/C with unconstrained marker covariance (`Mi.G.constraint == false`), JWAS calls `MTBayesABC_block!`.
- If `Mi.G.constraint == true`, JWAS uses `megaBayesABC!` (non-block path).
- In current multi-trait block mode, the update is sampler-I style (trait-wise `δ` updates from `BigPi[d0]`/`BigPi[d1]`) rather than the non-block sampler-I/II dispatcher.
- Current implementation note in source: this option is intended for one genotype category.
- In current implementation, numeric `fast_blocks` should satisfy `block_size < nMarkers` (chain-length scaling indexes the second block start).

## Single-Trait BayesC Without Blocks

In the standard BayesC update (`BayesABC!`), each marker `j` is updated one-by-one:

1. Compute conditional posterior terms (`rhs`, `lhs`, `gHat`) from current `yCorr`.
2. Compute marker inclusion probability `probDelta1`.
3. Sample `δ[j]` (include/exclude) and sample/update `β[j]`, `α[j]`.
4. Immediately update `yCorr` using marker column `x_j`.

So `yCorr` is updated every marker.

## Single-Trait Block BayesC in JWAS

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

- inner repeats are set to `nreps = block_size`;
- outer `chain_length` is reduced by approximately `block_size`.

This keeps effective marker-update work on a similar scale while moving much of the per-marker work from `nObs`-length operations to `block_size`-length operations.

### Detailed Comparison with BayesR3

Reference paper: BayesR3 (Communications Biology, 2022), DOI: `10.1038/s42003-022-03624-1`.

JWAS uses BayesC (not BayesR), but the block linear-algebra strategy closely follows the same blocked Gibbs pattern.

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
| Inner-repeat schedule | Uses nominal block size as fixed repeat count across blocks | `nreps = current block_size` | Last short block may receive fewer inner repeats |
| Outer-loop scheduling | Described as a fixed block sweep schedule | JWAS also rescales outer `chain_length` by block size | Compare effective updates, not just outer iterations |
| Scope | BayesR algorithm | JWAS block path is wired to BayesA/B/C marker samplers, with multi-trait block updates implemented in sampler-I style | Strategy reused in a different Bayesian alphabet member |

#### Scheduling detail (explicit)

JWAS sets `nreps` equal to the current block size.

- Full blocks: `nreps` equals the nominal block size.
- Final short block: `nreps` is smaller than full blocks.

In the BayesR3 description, `nreps` is treated as fixed by the nominal block size for all blocks, including the final short block.

This difference changes the number of within-block Gibbs sweeps for short blocks, but not the core residual/RHS block-update identity.

## Algorithm Comparison

| Aspect | Standard BayesC (`BayesABC!`) | Block BayesC (`BayesABC_block!`) | Practical effect |
| --- | --- | --- | --- |
| Update unit | One marker at a time | One block, then markers inside block | Better cache locality in block path |
| `yCorr` updates | Every marker | Once per block | Fewer full-length vector updates |
| Main per-marker linear algebra size | `nObs` | `block_size` (inside block cache) | Lower per-marker arithmetic cost |
| Extra precompute | Minimal | `X_b'R^{-1}X_b` for all blocks (RHS computed on demand) | More startup work |
| Extra memory | Minimal | Stores block Gram matrices (`XpRinvX`) and block workspaces | Higher memory than non-block, lower than older persisted-`XRinvArray` design |
| Chain behavior in current implementation | Direct `chain_length` | Inner repeats + outer chain scaling | Compare runs by effective updates, not only outer iterations |

## Computational Complexity

Use the notation:

- `N`: number of records (`nObs`)
- `P`: number of markers (`nMarkers`)
- `b`: nominal block size
- `B = ceil(P/b)`: number of blocks
- `L`: standard (non-block) chain length

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

JWAS rescales outer iterations to approximately `m = floor(L/b)`, so the main total cost is:

- `O((L/b) * (NP + P b^2)) = O(LP(N/b + b))`

### BayesR3 (block strategy and paper fit)

BayesR3 uses the same blocked-update strategy family (block RHS, in-block updates using a block Gram matrix, then a block-exit residual update), but it is a BayesR mixture model rather than BayesC. This changes constants (more mixture-state work per marker), not the core block linear-algebra pattern.

**Operation-count view (dense blocked implementation):**

- BayesR3 runs a *nominal* number of inner cycles `n` per block, and the paper recommends `n` be equal to the (nominal) block size `b`.
- With `n = b`, the leading operation-count terms match the same block strategy family as JWAS: total work scales like `O(LP(N/b + b))`.

**Paper runtime fit (Fig. 5):**

- The BayesR3 paper reports an empirical timing model where processing time per SNP is proportional to `(N + b)/b = N/b + 1`.
- This is a fit to measured runtime for their implementation/hardware and is not a formal asymptotic operation-count derivation. It effectively treats the in-block work (the `+b`-type term) as a small constant relative to the `N/b` term in that regime.

### Practical differences in complexity interpretation

- JWAS and BayesR3 share the same blocked-update strategy family, but are not identical in constants/scheduling details.
- JWAS uses `nreps = current block_size` for each block (including the final short block).
- BayesR3 describes using the nominal block repeat count for all blocks, including the final short block.

### Numerical example (`N=200,000`, `P=2,000,000`)

Assume `fast_blocks=true`, so JWAS uses `b = floor(sqrt(N)) = 447`.
Then:

- `B = ceil(P/b) = ceil(2,000,000/447) = 4,475` blocks
- Standard BayesC total scaling: `O(LPN) = O(L * 2,000,000 * 200,000)`
- JWAS block BayesC main scaling: `O(LP(N/b + b)) = O(L * 2,000,000 * (200,000/447 + 447))`
- BayesR3 paper timing fit: runtime per SNP `∝ (N + b)/b`, so total runtime `∝ L * 2,000,000 * (200,000/447 + 1)`

So the per-`LP` coefficients are:

- Standard BayesC: `200,000`
- JWAS block BayesC operation-count: `~894.4` (from `N/b + b`)
- BayesR3 paper fit: `~448.4` (from `N/b + 1`)

This implies:

- JWAS block vs standard: `~224x` lower
- BayesR3 paper fit vs standard: `~446x` lower
- The apparent `~2.0x` gap between `~894` and `~448` is not an apples-to-apples complexity comparison: it is the difference between an operation-count model (`N/b + b`) and an empirical runtime fit (`N/b + 1`).

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
