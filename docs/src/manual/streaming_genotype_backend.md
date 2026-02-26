# Streaming Genotype Backend (BayesC MVP)

This page documents the new opt-in streaming genotype backend for large-data
original BayesC in JWAS.

## Status and Design

- Dense loading remains the default and primary path.
  - `get_genotypes(...; storage=:dense)` (default)
  - Existing dense workflows are unchanged.
- Streaming loading is additive and opt-in.
  - `get_genotypes(...; storage=:stream)`
  - Current status: experimental MVP for original BayesC.

## Supported Scope (MVP)

- Single-trait analysis only
- Method: `BayesC` only
- `fast_blocks=false` only
- Unit residual weights only
- `double_precision=false` only
- Complete genomic data only (no single-step)
- Exact genotype/phenotype ID match and order required
- `outputEBV` and `output_heritability` are disabled in MVP streaming mode

## Public API

### `prepare_streaming_genotypes(...)`

One-time converter from text genotype file to packed marker-major 2-bit files.

```julia
prefix = prepare_streaming_genotypes(
    "genotypes.csv";
    output_prefix = "genotypes_stream",
    separator = ',',
    header = true,
    missing_value = 9.0,
    quality_control = true,
    MAF = 0.01,
    center = true,
    conversion_mode = :lowmem,   # :lowmem | :dense | :auto
    auto_dense_max_bytes = 2^30, # used only when conversion_mode=:auto
    tmpdir = nothing,          # optional temp workspace for conversion spool
    cleanup_temp = true,       # remove temporary conversion files after success
    disk_guard_ratio = 0.9,    # fail fast if required bytes exceed ratio * free bytes
)
```

Output artifacts with the chosen `prefix`:

- `prefix.jgb2`: packed genotype payload (2-bit marker-major)
- `prefix.meta`: metadata manifest
- `prefix.obsid.txt`: individual IDs
- `prefix.markerid.txt`: marker IDs (after QC)
- `prefix.mean.f32`: per-marker means
- `prefix.xpRinvx.f32`: per-marker `x'x` for unit-weight BayesC path
- `prefix.afreq.f32`: per-marker allele frequencies

Conversion notes:

- `conversion_mode=:lowmem` keeps conversion out-of-core (default).
- `conversion_mode=:dense` uses an in-memory dense conversion path (faster for small files, not RAM-safe for large files).
- `conversion_mode=:auto` chooses `:dense` when estimated dense bytes (`N*P*4`) fit under `auto_dense_max_bytes`; otherwise uses `:lowmem`.
- In low-memory mode, conversion performs a staged write path with temporary row-major spool + transpose.
- Temporary disk can approach one extra packed payload in low-memory mode; for very large `N, P`, place `tmpdir` on a high-capacity filesystem.

### How to use `conversion_mode`

Use one of these patterns:

```julia
# 1) Large files (safest): always low-memory conversion
prefix = prepare_streaming_genotypes("genotypes.csv";
                                     conversion_mode=:lowmem,
                                     tmpdir="/scratch/jwas_tmp")

# 2) Small files (fast path): force dense in-memory conversion
prefix = prepare_streaming_genotypes("genotypes.csv";
                                     conversion_mode=:dense)

# 3) Hybrid default: auto-select dense for small jobs, lowmem otherwise
prefix = prepare_streaming_genotypes("genotypes.csv";
                                     conversion_mode=:auto,
                                     auto_dense_max_bytes=2^30) # 1 GiB threshold
```

Practical rule: start with `conversion_mode=:auto`; use `:lowmem` explicitly for very large jobs or constrained RAM environments.

### `get_genotypes(...; storage=:stream)`

Load packed genotype backend without materializing dense `N x P` genotype matrix.

```julia
geno = get_genotypes(
    prefix,            # prefix or .meta/.jgb2 path
    1.0;
    method = "BayesC",
    estimatePi = true,
    storage = :stream,
)
```

Dense behavior is unchanged:

```julia
geno = get_genotypes("genotypes.csv", 1.0; method="BayesC")  # storage=:dense default
```

## End-to-End Example

```julia
using JWAS, CSV, DataFrames

phenotypes = CSV.read("phenotypes.csv", DataFrame)

# one-time conversion
prefix = prepare_streaming_genotypes("genotypes.csv";
                                     separator=',',
                                     header=true,
                                     quality_control=true,
                                     center=true)

# streaming load
global geno = get_genotypes(prefix, 1.0;
                            method="BayesC",
                            estimatePi=true,
                            storage=:stream)

model = build_model("trait1 = intercept + geno", 1.0)

out = runMCMC(model, phenotypes;
              chain_length=500,
              burnin=100,
              output_samples_frequency=50,
              output_folder="results_stream",
              outputEBV=false,
              output_heritability=false,
              seed=314,
              memory_guard=:off)
```

## Benchmark Snapshot

All numbers below were measured on:

- Julia `1.11.7`
- Apple M1 (8 CPU threads available)
- macOS arm64

### 1) Correctness: `simulated_omics` dense vs stream

Dataset:

- `/src/4.Datasets/data/simulated_omics`
- 3534 individuals
- 1000 SNPs (927 markers after QC)

Same seed/settings, single-trait BayesC:

| Metric | Value |
| --- | ---: |
| Marker effects correlation | `0.999999999978` |
| Marker max abs diff | `1.10e-6` |
| Marker mean abs diff | `1.65e-7` |
| Residual variance abs diff | `4.29e-7` |

Interpretation:

- Streaming and dense are numerically equivalent up to floating-point noise.

### 2) Large synthetic benchmark (`N=10,000`, `P=5,000`)

Synthetic CSV genotype file size: `~100 MB`.

`runMCMC` benchmark (`chain_length=25`, `burnin=5`):

| Mode | Real Time | Max RSS | Peak footprint |
| --- | ---: | ---: | ---: |
| Dense (load + run) | `24.19s` | `1.396 GB` | `1.352 GB` |
| Stream run (after prepare) | `20.70s` | `0.931 GB` | `0.727 GB` |
| Stream prepare (one-time) | `11.99s` | `1.367 GB` | `1.144 GB` |

Interpretation:

- Streaming reduced run-time memory significantly in this benchmark.
- Conversion has one-time cost and memory use.
- Absolute times are hardware- and IO-dependent.

### 2b) Conversion-phase benchmark (`N=10,000`, `P=5,000`)

Measured with `/usr/bin/time -l`:

| Step | Real Time | Max RSS |
| --- | ---: | ---: |
| `prepare_streaming_genotypes` | `~11.8s` | `~1.14 GB` |

Backend sizes:

| Artifact | Size |
| --- | ---: |
| packed payload (`.jgb2`) | `~12.5 MB` |
| temporary row-major spool (during conversion) | `~12.5 MB` |

Interpretation:

- Converter RAM stays bounded for this size and no dense matrix is materialized.
- Conversion needs temporary disk in addition to final payload.

### 3) Target-scale feasibility estimator (`N=500,000`, `P=2,000,000`)

Using `estimate_marker_memory(...)`:

| Path | Estimated marker-path memory |
| --- | ---: |
| Dense (`storage=:dense`, Float32) | `3.64 TiB` |
| Stream (`storage=:stream`, Float32) | `17.29 MiB` |

Interpretation:

- Dense original BayesC is generally infeasible at this scale on normal nodes.
- Streaming is designed to keep marker working memory near `O(N+P)`.

## Notes

- Streaming support for non-unit weights and `fast_blocks` is planned as follow-up work.
- For algorithmic memory/speed derivations, see
  [Handling Large Genotype Data Without Loading the Full Matrix into Memory](large_genotype_data_streaming.md).
