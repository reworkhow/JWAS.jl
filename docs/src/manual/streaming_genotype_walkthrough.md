# Streaming Genotype Workflow: A Conceptual Walkthrough

This page walks through the streaming genotype backend from start to finish, without code execution.
For API details and constraints, see [Streaming Genotype Backend](streaming_genotype_backend.md)
and [Handling Large Genotype Data](large_genotype_data_streaming.md).

## 1. The Problem

At very large scale (e.g., N=500,000 individuals × P=2,000,000 markers), a dense genotype matrix
in memory requires terabytes of RAM (~4 TB for Float32). The streaming backend avoids ever materializing
that matrix by keeping genotypes on disk in a compact format and decoding only what is needed,
one marker at a time.

## 2. Input: Text Genotype File

You start with a text genotype file (e.g., CSV):

- **First column**: individual IDs
- **Remaining columns**: marker genotypes (0, 1, 2)
- **Missing values**: coded as 9.0 (or your chosen value)

Example format:

```
ID,m1,m2,m3,...,m100
ind_001,0,1,2,...,1
ind_002,1,9,0,...,2   ← 9 = missing
...
```

## 3. One-Time Conversion: `prepare_streaming_genotypes`

This step reads the text file and writes a packed backend. It runs once, before any MCMC.

**What it does:**

1. **Read and process** the text file. The current implementation loads the full genotype matrix into memory during this step, so `prepare_streaming_genotypes` must run on a machine with enough RAM to hold the data. The memory savings apply to the MCMC phase, not to the conversion step.
2. **Quality control**: replace missing with column means; drop fixed loci and markers below the MAF threshold.
3. **Center** each marker column (optional).
4. **Precompute** per-marker x′x (for unit-weight BayesC).
5. **Pack** genotypes into 2-bit codes (0→00, 1→01, 2→10, missing→11) and write marker-major to disk.

**Output files:**

- `prefix.jgb2` — binary file containing 2-bit packed genotypes (0, 1, 2, missing) in marker-major layout; four individuals per byte
- `prefix.meta` — manifest (nObs, nMarkers, paths, etc.)
- `prefix.obsid.txt`, `prefix.markerid.txt` — individual and marker IDs
- `prefix.mean.f32`, `prefix.xpRinvx.f32`, `prefix.afreq.f32` — precomputed values

After QC, markers are stored in column order. For 50 individuals, each marker occupies ~13 bytes
(4 bits per individual, rounded up).

## 4. Loading: `get_genotypes(prefix; storage=:stream)`

Instead of loading a dense N×P matrix, this step:

1. **Reads** the manifest (the `prefix.meta` file: a small text file listing paths to the packed data and metadata such as nObs, nMarkers, stride_bytes, centered) and validates file sizes.
2. **Opens** the `.jgb2` file for reading.
3. **Loads** metadata into memory: means, xpRinvx, allele frequencies, IDs.
4. **Creates** a `Genotypes` struct with:
   - `genotypes` = empty matrix (no dense X)
   - `stream_backend` = handle to the packed file and metadata
   - `storage_mode` = `:stream`

## 5. MCMC: Per Iteration, Per Marker

For each MCMC iteration, BayesC updates markers one at a time.

**For each marker j:**

1. **Decode** marker j into a temporary buffer:
   - Seek to the marker's byte offset in the packed file.
   - Read one packed row.
   - Decode 2-bit codes → 0, 1, 2, or mean-imputed value.
   - Apply centering if the backend was built with `center=true`.
   - Write into a length-N buffer.

2. **Update** Gibbs:
   - Use the same formulas as dense mode.
   - Use `dot(buffer, yCorr)` for x′y.
   - Use precomputed `xpRinvx[j]` for x′x.
   - Sample δ, β, α.
   - Update `yCorr` with `yCorr += (old_α - new_α) * buffer`.

3. **Discard** the buffer; it is reused for the next marker.

So the working set is O(N) per marker (one buffer plus yCorr, α, β, δ), not O(N×P).

## 6. Initial yCorr

The full computation `yCorr = y - X*sol - X*α` is done **once, before the MCMC loop** (using starting values for sol and α). After that, yCorr is maintained by incremental updates: each iteration adds/subtracts X*sol for fixed-effect updates, and each marker update applies `yCorr += (old_α - new_α) * x_j` in place.

Code: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`, lines 108–124 (inside `MCMC_BayesianAlphabet`, before the iteration loop).

**Dense:** `yCorr = y - X*sol - X*α` (matrix multiply).

**Streaming:** `streaming_mul_alpha!` replaces `X*α`:
- Loop over markers j.
- For each j with α[j] ≠ 0: decode marker j into a buffer, add `α[j] * buffer` into the output vector.
- No dense X is ever formed.

## 7. MVP Constraints

- Single-trait analysis only
- Method: BayesC only
- `fast_blocks=false` only
- Unit residual weights only
- Float32 only (`double_precision=false`)
- Complete genomic data only (no single-step)
- Phenotype IDs must match genotype IDs exactly and in the same order (no alignment)
- `outputEBV` and `output_heritability` are disabled

## 8. End-to-End Flow

```
Text CSV → prepare_streaming_genotypes() → .jgb2 + metadata
                    ↓
        get_genotypes(prefix; storage=:stream)
                    ↓
        Genotypes struct (no dense X)
                    ↓
        runMCMC()
                    ↓
Per iteration:
  Per marker j:
    decode_marker!(buffer, backend, j) → buffer of length N
    bayesabc_update_marker!(buffer, yCorr, α, β, δ, ...)
```

## 9. Memory Comparison

| Component        | Dense (N×P) | Streaming |
|-----------------|-------------|-----------|
| Genotype matrix  | N×P×4 bytes | 0         |
| Per-marker buffer | 0        | N×4 bytes |
| xpRinvx         | P×4 bytes  | P×4 bytes |
| **Total**       | ~4×N×P     | ~O(N + P) |

At large N and P, streaming keeps memory roughly constant in N×P, at the cost of decode and I/O per marker per iteration.
