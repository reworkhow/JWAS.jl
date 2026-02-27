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

## Cluster BayesC Benchmarks (Fast-Block vs Non-Block)

This folder also includes Slurm-oriented scripts used for the benchmark page:

- `benchmarks/jwas_full_benchmark.jl`
- `benchmarks/jwas_full_benchmark.sbatch`
- `benchmarks/jwas_full_benchmark_v2.jl`
- `benchmarks/jwas_full_benchmark_v2.sbatch`
- `benchmarks/jwas_nonblock_benchmark.jl`
- `benchmarks/jwas_nonblock_benchmark.sbatch`

### What each script does

- `jwas_full_benchmark.*`: initial fast-block benchmark run.
- `jwas_full_benchmark_v2.*`: fast-block benchmark with warmup + replicates (recommended).
- `jwas_nonblock_benchmark.*`: standard BayesC (`fast_blocks=false`) benchmark with warmup + replicates.

### Before submitting

These `.sbatch` files are templates from a real cluster run and contain cluster-specific paths.
Update the following in each `.sbatch` file:

- partition/account (`#SBATCH -p`, `#SBATCH -A`)
- output/error paths (`#SBATCH -o`, `#SBATCH -e`)
- Julia module line (if your cluster differs)
- `--project=...` path in the final `julia` command
- benchmark script path in the final `julia` command
- output text path env var (`BENCH_TXT`)

### Example usage

From repository root:

```bash
# 1) Fast-block (recommended replicated run)
sbatch benchmarks/jwas_full_benchmark_v2.sbatch

# 2) Non-block benchmark
sbatch benchmarks/jwas_nonblock_benchmark.sbatch
```

Monitor:

```bash
squeue -u $USER
tail -f /path/to/slurm_output.out
```

Find final estimate:

```bash
rg "SUMMARY_EST|EXTRAPOLATED|RESULTS_WRITTEN" /path/to/slurm_output.out
cat /path/to/benchmark_results.txt
```

### Runtime controls

These benchmarks are controlled via environment variables in the `.sbatch` files:

- `BENCH_N`: number of individuals
- `BENCH_P_LIST`: calibration marker counts (comma-separated)
- `BENCH_L_LIST`: calibration chain lengths (comma-separated)
- `BENCH_REPS`: replicate count (in v2/non-block scripts)
- `WARMUP_N`, `WARMUP_P`, `WARMUP_L`: warmup size
- `TARGET_P`, `TARGET_L`: extrapolation target

If you change target scale, keep calibration sizes feasible for your node memory and wall-time limits.

### Farm preset (tested on 2026-02-26)

The settings below were used for the published `Benchmark` doc page and can be used as a known-good starting point on `farm`.

Common Slurm settings:

- partition: `low`
- account: `publicgrp`
- nodes/tasks: `--nodes=1 --ntasks=1`
- CPU/memory: `--cpus-per-task=64 --mem=220G`
- Julia module: `julia/1.12.1`
- project path: `/home/qtlcheng/Github/JWAS.jl`

Fast-block (replicated, recommended):

- script: `benchmarks/jwas_full_benchmark_v2.sbatch`
- walltime: `04:00:00`
- env:
  - `BENCH_N=50000`
  - `BENCH_P_LIST=100000,200000`
  - `BENCH_L_LIST=223,446`
  - `BENCH_REPS=2`
  - `WARMUP_N=2000`, `WARMUP_P=5000`, `WARMUP_L=50`
  - `TARGET_P=2000000`, `TARGET_L=2000`
- completed benchmark job: `30951323` (node `cpu-10-63`)

Non-block (standard BayesC):

- script: `benchmarks/jwas_nonblock_benchmark.sbatch`
- walltime: `06:00:00`
- env:
  - `BENCH_N=50000`
  - `BENCH_P_LIST=100000,200000`
  - `BENCH_L_LIST=2,4`
  - `BENCH_REPS=2`
  - `WARMUP_N=2000`, `WARMUP_P=2000`, `WARMUP_L=5`
  - `TARGET_P=2000000`, `TARGET_L=2000`
- completed benchmark job: `30951596` (node `cpu-6-91`)
