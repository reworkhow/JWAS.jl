# Block BayesC (`fast_blocks`)

This page explains how JWAS implements block updates for BayesC marker sampling and how the block path changes speed and memory usage.

## When This Path Is Used

Block updates are enabled with:

```julia
out = runMCMC(model, phenotypes; fast_blocks=true)
# or provide a fixed block size
out = runMCMC(model, phenotypes; fast_blocks=64)
```

- If `fast_blocks=true`, JWAS chooses `block_size = floor(sqrt(nObs))`.
- In single-trait BayesA/B/C, JWAS calls `BayesABC_block!`.
- In multi-trait BayesA/B/C, JWAS calls `MTBayesABC_block!`.
- Current implementation note in source: this option is intended for one genotype category.

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
- `X_b'R^{-1}` (`XRinvArray`)
- `X_b'R^{-1}X_b` (`XpRinvX`, block Gram matrix)

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

## Algorithm Comparison

| Aspect | Standard BayesC (`BayesABC!`) | Block BayesC (`BayesABC_block!`) | Practical effect |
| --- | --- | --- | --- |
| Update unit | One marker at a time | One block, then markers inside block | Better cache locality in block path |
| `yCorr` updates | Every marker | Once per block | Fewer full-length vector updates |
| Main per-marker linear algebra size | `nObs` | `block_size` (inside block cache) | Lower per-marker arithmetic cost |
| Extra precompute | Minimal | `X_b'R^{-1}` and `X_b'R^{-1}X_b` for all blocks | More startup work |
| Extra memory | Minimal | Stores block matrices | Higher memory footprint |
| Chain behavior in current implementation | Direct `chain_length` | Inner repeats + outer chain scaling | Compare runs by effective updates, not only outer iterations |

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
| Extra `XRinvArray` in block mode | `~953.7 MiB` |
| Extra `XpRinvX` in block mode | `~13.4 MiB` |

Interpretation:

- Block mode can be much faster for large `nObs`, because heavy per-marker operations are shifted to block-level cached operations.
- Block mode uses more memory, mainly from `XRinvArray`.
- Practical speedup is typically below theoretical arithmetic speedup due to random branching, allocation overhead, and BLAS/runtime effects.

## Practical Guidance

1. Start with `fast_blocks=true` for large marker sets and enough RAM.
2. If memory is tight, set a smaller numeric block size (e.g., `32` or `64`) and benchmark.
3. If speed gain is small, try a few block sizes and choose based on wall time + memory headroom.

