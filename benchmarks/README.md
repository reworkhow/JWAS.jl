# Benchmarks

This folder contains reproducible benchmark scripts for large-data streaming
BayesC performance and memory studies.

## Conversion Benchmark Notes

`prepare_streaming_genotypes(...)` is a one-time out-of-core conversion stage.
It does not load dense `N x P` into RAM, but it uses temporary disk during
row-major spool + transpose.

Observed (`N=10,000`, `P=5,000`, `/usr/bin/time -l`):

- conversion wall time: `~11.8s`
- conversion max RSS: `~1.14 GB`
- final packed payload (`.jgb2`): `~12.5 MB`
- temporary spool during conversion: `~12.5 MB`

## Streaming Large Benchmark

Script:

- `benchmarks/streaming_large_benchmark.jl`

What it does:

1. Optionally creates a synthetic sparse streaming backend at target `N` and `P`.
2. Loads backend through JWAS `storage=:stream`.
3. Benchmarks marker-update-like loops for configurable sweep counts.
4. Samples process RSS and writes CSV report.
5. Reports dense vs stream memory estimates from `estimate_marker_memory`.

### Quick sanity run

```bash
julia --project benchmarks/streaming_large_benchmark.jl \
  --create-backend=true \
  --force-create=true \
  --n-obs=10000 \
  --n-markers=5000 \
  --prefix=/tmp/jwas_stream_10k_5k \
  --sweeps=1,5 \
  --sample-markers-per-sweep=1000 \
  --csv-out=/tmp/stream_bench_10k_5k.csv
```

### Target-scale benchmark (`N=500,000`, `P=2,000,000`)

Full-sweep mode:

```bash
julia --project benchmarks/streaming_large_benchmark.jl \
  --create-backend=true \
  --n-obs=500000 \
  --n-markers=2000000 \
  --prefix=/scratch/jwas_stream_500k_2m \
  --sweeps=1,5,10 \
  --sample-markers-per-sweep=0 \
  --csv-out=/scratch/jwas_stream_500k_2m_fullsweep.csv
```

Sampled mode (faster local projection run):

```bash
julia --project benchmarks/streaming_large_benchmark.jl \
  --create-backend=true \
  --n-obs=500000 \
  --n-markers=2000000 \
  --prefix=/scratch/jwas_stream_500k_2m \
  --sweeps=1,5,10 \
  --sample-markers-per-sweep=20000 \
  --csv-out=/scratch/jwas_stream_500k_2m_sampled.csv
```

### Key arguments

- `--create-backend=true|false`: create synthetic backend if needed.
- `--force-create=true|false`: overwrite/recreate backend files.
- `--n-obs`, `--n-markers`: backend dimensions.
- `--prefix`: backend file prefix.
- `--sweeps=1,5,10`: sweep counts to benchmark.
- `--sample-markers-per-sweep=K`: `0` means full marker sweep; `K>0` means sampled subset per sweep with projection.
- `--rss-sample-interval-sec`: RSS sampling interval.
- `--csv-out`: CSV output path.
