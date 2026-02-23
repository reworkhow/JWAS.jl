# Memory Usage (BayesA/B/C Marker Paths)

This page summarizes memory usage for JWAS marker samplers, with emphasis on the non-block path (`fast_blocks=false`) and how block mode changes memory requirements.

## Notation

- `N`: number of records (`nObs`)
- `P`: number of markers (`nMarkers`)
- `b`: nominal block size
- `B = ceil(P / b)`: number of blocks
- `s_i`: size of block `i`, with `sum_i s_i = P`
- `t`: bytes per stored value (`t=4` for `Float32`, `t=8` for `Float64`)

Main allocations are created in `GibbsMats(...)` (`src/1.JWAS/src/markers/tools4genotypes.jl`).

## Non-Block Mode (`fast_blocks=false`)

### What is stored

- `X` (genotype matrix): `N x P`
- `xArray`: column views of `X` (metadata only, no data copy)
- `xpRinvx`: length `P`
- `xRinvArray`:
  - aliases `xArray` when `Rinv == ones(length(Rinv))` (no extra `N x P` copy)
  - materializes `[x .* Rinv for x in xArray]` when weights are non-unit (extra `N x P`)

There are also marker-state vectors (`alpha`, `beta`, `delta`, posterior means, etc.) that scale as `O(P)` and are usually much smaller than `N x P` matrices.

### Approximate memory formulas

- Unit weights (default when no `weights` column is used):  
`Mem_nonblock_unit ~= t * (N*P + P) + O(P*t)`
- Non-unit weights:  
`Mem_nonblock_nonunit ~= t * (2*N*P + P) + O(P*t)`

Interpretation: non-block memory is dominated by one copy of `X` (unit weights) or two copies (non-unit weights).

## Block Mode Additions (`fast_blocks != false`)

Block mode keeps all non-block structures and adds:

- `XArray`: block views of `X` (metadata only)
- `XRinvArray`: per-block dense matrices `X_b' * Diagonal(Rinv)`; total elements `N*P`
- `XpRinvX`: per-block Gram matrices; total elements `sum_i s_i^2` (approximately `P*b` for near-uniform blocks)

Approximate totals:

- Unit weights:  
`Mem_block_unit ~= t * (2*N*P + sum_i s_i^2 + P) + O(P*t)`
- Non-unit weights:  
`Mem_block_nonunit ~= t * (3*N*P + sum_i s_i^2 + P) + O(P*t)`

## Worked Example (`N=500,000`, `P=2,000,000`)

Assume `fast_blocks=true`, so `b=floor(sqrt(N))=707`.

- `B = ceil(P/b) = 2,829`
- `sum_i s_i^2 = 1,413,937,788`

### Original non-block version (`fast_blocks=false`)

This is the base/original sampler memory footprint without block matrices.

| Component | Elements | Float32 | Float64 |
| --- | ---: | ---: | ---: |
| `X` | `N*P = 1,000,000,000,000` | `4.00 TB` (`3.64 TiB`) | `8.00 TB` (`7.28 TiB`) |
| `xRinvArray` (unit weights) | alias of `xArray` (no data copy) | `0` | `0` |
| `xRinvArray` (non-unit weights) | `N*P` | `4.00 TB` (`3.64 TiB`) | `8.00 TB` (`7.28 TiB`) |
| `xpRinvx` | `P` | `8.0 MB` (`7.63 MiB`) | `16.0 MB` (`15.26 MiB`) |

Non-block totals:

| Non-block case | Float32 | Float64 |
| --- | ---: | ---: |
| Unit weights | `~4.00 TB` | `~8.00 TB` |
| Non-unit weights | `~8.00 TB` | `~16.00 TB` |

### Extra objects introduced by block mode (`fast_blocks != false`)

Component sizes:


| Component                            | Elements                  | Float32                | Float64                  |
| ------------------------------------ | ------------------------- | ---------------------- | ------------------------ |
| `X`                                  | `N*P = 1,000,000,000,000` | `4.00 TB` (`3.64 TiB`) | `8.00 TB` (`7.28 TiB`)   |
| `xRinvArray` (only non-unit weights) | `N*P`                     | `4.00 TB` (`3.64 TiB`) | `8.00 TB` (`7.28 TiB`)   |
| `XRinvArray` (block mode)            | `N*P`                     | `4.00 TB` (`3.64 TiB`) | `8.00 TB` (`7.28 TiB`)   |
| `XpRinvX` (block mode)               | `sum_i s_i^2`             | `5.66 GB` (`5.27 GiB`) | `11.31 GB` (`10.53 GiB`) |
| `xpRinvx`                            | `P`                       | `8.0 MB` (`7.63 MiB`)  | `16.0 MB` (`15.26 MiB`)  |


Approximate totals (dominated by `N*P` terms):


| Mode                        | Float32     | Float64     |
| --------------------------- | ----------- | ----------- |
| Non-block, unit weights     | `~4.00 TB`  | `~8.00 TB`  |
| Non-block, non-unit weights | `~8.00 TB`  | `~16.00 TB` |
| Block, unit weights         | `~8.01 TB`  | `~16.01 TB` |
| Block, non-unit weights     | `~12.01 TB` | `~24.01 TB` |


## Practical Takeaways

1. At very large `N` and `P`, dense in-memory genotype matrices dominate memory in both non-block and block modes.
2. Non-block mode is memory-cheapest when weights are unit (`xRinvArray` aliases `xArray`).
3. Block mode adds one extra `N x P`-scale matrix (`XRinvArray`) plus a smaller `XpRinvX` term.
4. For datasets like `N=500k, P=2M`, dense storage is multi-terabyte and typically requires a different data representation strategy.
