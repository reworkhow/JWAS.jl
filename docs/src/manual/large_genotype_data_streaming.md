# Handling Large Genotype Data Without Loading the Full Matrix into Memory

This page summarizes scalability limits and representation options for the standard
original BayesC path (`fast_blocks=false`), with focus on very large datasets such as
`N=500,000` individuals and `P=2,000,000` markers.

## Scope

- Model path: standard BayesC, non-block marker updates.
- Objective: understand memory and speed constraints, and compare storage approaches.
- This page does not change BayesC posterior equations; it compares data representation and runtime behavior.

## Current JWAS Status

- Dense loading remains the default and primary path:
  - `get_genotypes(...; storage=:dense)` (default)
  - Existing dense workflows are unchanged.
- Streaming loading is additive and opt-in:
  - `get_genotypes(...; storage=:stream)`
  - Status: experimental MVP for large-data BayesC.
- Streaming conversion is out-of-core:
  - `prepare_streaming_genotypes(...)` no longer materializes dense `N x P` in RAM.
  - Conversion uses temporary disk plus final packed output (set `tmpdir` for placement).

### Streaming MVP workflow

```julia
# 1) one-time conversion to packed backend
prefix = prepare_streaming_genotypes(\"genotypes.csv\";
                                     separator=',',
                                     header=true,
                                     quality_control=true,
                                     center=true,
                                     tmpdir=nothing,
                                     cleanup_temp=true,
                                     disk_guard_ratio=0.9)

# 2) load packed backend (no dense N x P matrix in memory)
geno = get_genotypes(prefix, 1.0;
                     method=\"BayesC\",
                     storage=:stream)
```

### Streaming MVP constraints

- single-trait only
- method = `BayesC` only
- `fast_blocks=false` only
- unit residual weights only
- `double_precision=false` only
- complete genomic data only (no single-step)
- exact genotype/phenotype ID match and order required
- `outputEBV`/`output_heritability` are disabled in MVP streaming mode

## Notation

- `N`: number of records (`nObs`)
- `P`: number of markers (`nMarkers`)
- `t`: bytes per stored value (`4` for `Float32`, `8` for `Float64`)
- `L`: chain length (effective non-block iterations)
- `c`: number of markers decoded/loaded per chunk in an out-of-core path
- `B_io`: effective sustained storage throughput (GB/s)

## Baseline Complexity (Original BayesC)

For single-trait original BayesC, each MCMC iteration updates markers one-by-one.
The dominant per-iteration cost scales as:

- `O(NP)`

Total over `L` iterations:

- `O(LNP)`

At very large `N` and `P`, runtime is dominated by repeated full marker sweeps.

### Rough operation count

For each marker update in the single-trait loop:

- one `dot(x_j, yCorr)` term (`~2N` floating-point ops),
- one residual update `yCorr += (old-new)*x_j` (`~2N` ops),
- small `O(1)` scalar work.

A practical back-of-the-envelope is:

- per iteration: `~4NP` flops
- total: `~4LNP` flops

For `N=500,000`, `P=2,000,000`:

- `~4e12` flops per iteration (before overhead/branching/cache effects)

## Baseline Memory (Original BayesC)

Main memory terms for one genotype category:

- Dense genotype matrix `X`: `N x P`
- `xpRinvx`: length `P`
- `xRinvArray`:
  - unit weights: aliases `xArray` (no extra `N x P` data copy)
  - non-unit weights: extra `N x P` copy

Approximate totals:

- Unit weights: `Mem_nonblock_unit ~= t * (N*P + P)`
- Non-unit weights: `Mem_nonblock_nonunit ~= t * (2*N*P + P)`

In practice, `N*P` dominates.

## Worked Example (`N=500,000`, `P=2,000,000`)

- `N*P = 1,000,000,000,000`
- 2-bit packed genotype payload size (for comparison): `N*P/4 = 250,000,000,000` bytes

### Memory totals for original BayesC

| Case | Float32 | Float64 |
| --- | ---: | ---: |
| Unit weights | `~4.00 TB` (`3.64 TiB`) | `~8.00 TB` (`7.28 TiB`) |
| Non-unit weights | `~8.00 TB` (`7.28 TiB`) | `~16.00 TB` (`14.55 TiB`) |

## Out-of-Core Working-Set Math (Original BayesC)

If genotypes are streamed by chunks of `c` markers and decoded into Float32:

- chunk buffer: `N*c*4` bytes
- marker-state vectors (`α`, `β`, `δ`): `O(P)` (small relative to chunk for large `N`)
- response/residual vectors (`y`, `yCorr`): `O(N)` (also small relative to chunk)

Approximate runtime working set:

- `Mem_working ~= 4*N*c + O(P) + O(N)` bytes

Example at `N=500,000`:

| Chunk size `c` | `X_chunk` buffer (Float32) |
| ---: | ---: |
| 128 | `~256 MB` |
| 256 | `~512 MB` |
| 512 | `~1.02 GB` |
| 1024 | `~2.05 GB` |

This is why out-of-core design can be RAM-feasible even when dense in-memory is not.

## Representation Approaches for Original BayesC

### 1. Dense In-Memory (current baseline)

- Store dense `X` in RAM.
- Best arithmetic locality.
- Requires multi-terabyte RAM at this scale.
- No extra decode overhead.

### 2. Dense Out-of-Core (`mmap`-style dense backend)

- Dense storage on disk, streamed by marker sweep.
- Reduces RAM requirement, but disk remains dense-scale.

For the example:
- Dense file size (`Float32`): `~4.00 TB`

Pros:
- lowest engineering risk
- minimal change to current math path

Cons:
- very large disk footprint
- heavy I/O per iteration
- usually bottlenecked by storage bandwidth, not compute

### 3. Native Bit-Packed Backend (recommended long-term)

- Store genotypes in compact 2-bit representation.
- Decode on demand into small work buffers during marker updates.

For the example:
- Packed payload: `~250 GB` (plus metadata/index overhead)

Pros:
- major reduction in disk footprint and read traffic
- scalable foundation for large `N, P`
- can target predictable RAM via chunk size `c`

Cons:
- higher engineering effort (decoder, metadata, validation)
- decode overhead during runtime

### 4. External Adapter Backend (e.g., PLINK-style reader)

- Reuse existing external compact format and stream into BayesC update loops.

Pros:
- interoperability with existing data pipelines

Cons:
- parser/dependency complexity
- performance depends on adapter implementation and access pattern

## Fast-Block Implementation (`fast_blocks`) in Large-Data Context

Although this page focuses on original BayesC (`fast_blocks=false`), large-data
planning usually compares it to `fast_blocks=true`.

Let:

- `b`: nominal block size
- `s_i`: block sizes with `sum_i s_i = P`
- `S = sum_i s_i^2` (near-uniform approximation: `S ~ P*b`)

Current block path (after removing persistent `XRinvArray`) is approximately:

- Unit weights:
  - `Mem_block_unit ~= t * (N*P + S + P) + O(b*t)`
- Non-unit weights:
  - `Mem_block_nonunit ~= t * (2*N*P + S + P) + O(b*t)`

Computation scales approximately as:

- per outer iteration: `O(NP + sum_i s_i^3)` (`~O(NP + P*b^2)` for near-uniform blocks)
- with JWAS outer-loop rescaling (`m ~ L/b`): `O(LP(N/b + b))`

So `fast_blocks` can substantially reduce arithmetic relative to original BayesC,
while adding mainly the `S` term for `XpRinvX`.

## What To Watch for `fast_blocks`

1. Block-size tradeoff:
   - smaller `b`: less memory, weaker speedup
   - larger `b`: higher `XpRinvX` memory and heavier precompute
2. Effective iteration accounting:
   - compare runs by effective updates, not only reported outer iterations
3. Final short block:
   - last block gets fewer inner repeats because `nreps = block_size` per block
4. Startup/precompute cost:
   - `XpRinvX` build can dominate startup on very large `N,P`
5. Multi-trait specifics:
   - current multi-trait block mode behavior differs from the full non-block sampler dispatcher
6. Numerical reproducibility:
   - equivalent algebraic refactors may change floating-point roundoff (especially `Float32`)
7. Weighting caveat:
   - block `XRinvArray` is no longer persisted, but non-unit weighted non-block `xRinvArray` remains a separate optimization topic

## Speed and I/O Considerations

Original BayesC still has `O(LNP)` compute scaling regardless of representation.
Representation changes mostly affect:

- feasible RAM footprint
- disk footprint
- I/O volume per marker sweep

For out-of-core paths, each iteration needs at least one full marker sweep.
Thus, reducing per-sweep data size (e.g., dense vs 2-bit packed) directly reduces I/O pressure.

### I/O lower-bound model

Ignoring compute and decode for a lower bound:

- `T_io_per_iter >= Bytes_per_sweep / B_io`
- `T_io_total >= L * Bytes_per_sweep / B_io`

For the worked example (`Float32` dense vs 2-bit packed):

- dense sweep bytes: `~4,000 GB`
- packed sweep bytes: `~250 GB`

If `B_io = 3 GB/s`:

| Storage form | I/O lower bound per iteration |
| --- | ---: |
| Dense Float32 | `~22.2 min` |
| 2-bit packed | `~1.39 min` |

At `L=1000` iterations (I/O lower bound only):

| Storage form | I/O lower bound total |
| --- | ---: |
| Dense Float32 | `~15.4 days` |
| 2-bit packed | `~23.1 hours` |

These are lower bounds; actual wall time is higher once decode and BayesC compute are included.

## Side-by-Side Summary

| Approach | RAM feasibility at 500k x 2M | Disk footprint | I/O pressure | Engineering effort |
| --- | --- | --- | --- | --- |
| Dense in-memory | very poor | low (if already in RAM source) | low | low |
| Dense out-of-core (`mmap`) | good | very high (`~4 TB`, Float32) | very high | low |
| Native 2-bit packed | good | moderate (`~250 GB` payload) | much lower | high |
| External compact adapter | good | format-dependent (often compact) | lower than dense | medium-high |

## Validation Details Needed for Large-Data Streaming

For any out-of-core backend, validate:

1. numerical equivalence vs dense baseline on small/medium data (`α`, `δ`, `yCorr` traces),
2. missing/imputation and centering consistency,
3. deterministic behavior with fixed RNG seed,
4. chunk-size invariance (same results for different `c` within FP tolerance),
5. sustained throughput and memory caps under long chains.

## Practical Takeaways

1. At `N=500k, P=2M`, original BayesC dense in RAM is not practical on typical hardware.
2. Out-of-core dense storage improves RAM feasibility but remains very heavy in disk/I/O.
3. A compact backend (native bit-packed or external compact format) is the practical direction for large-scale original BayesC.
4. Even with compact storage, runtime remains fundamentally `O(LNP)`; representation helps feasibility, not asymptotic compute.
