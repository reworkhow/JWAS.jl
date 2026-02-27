# Benchmark

This page records a real benchmark run comparing `runMCMC(...; fast_blocks=true)` and `runMCMC(...; fast_blocks=false)` for the same target-scale question:

- `N = 50,000` individuals
- target `P = 2,000,000` markers
- target chain length `L = 2000`

Date of run: February 26, 2026.

## Goal

Estimate wall-clock MCMC time for a single-trait BayesC run at target scale, and compare fast-block vs standard BayesC.

## Environment

- Platform: UCD farm cluster (Slurm)
- Julia: `1.12.1`
- Model: `y1 = intercept + geno`
- Data: synthetic dense `Float32` genotype matrix and synthetic phenotype
- Settings:
  - `method="BayesC"`
  - `memory_guard=:warn`
  - `outputEBV=false`
  - `output_heritability=false`

## Fast-Block Benchmark (`fast_blocks=true`)

### Design

- Slurm job: `30951323`
- Warmup run before timed cases
- Timed grid:
  - `P ∈ {100,000, 200,000}`
  - input `chain_length ∈ {223, 446}`
  - `REPS = 2`
- Because `fast_blocks=true` and `N=50,000`:
  - `block_size = floor(sqrt(50000)) = 223`
  - outer iterations are `floor(chain_length / 223)` (1 or 2 in timed runs)
- Extrapolation:
  1. For each `P`, fit `t_mcmc = a + b * outer`
  2. Convert target to outer iterations: `floor(2000 / 223) = 8`
  3. Fit target time vs `P` using `P = 100,000` and `P = 200,000`, then extrapolate to `P = 2,000,000`

### Measured `t_mcmc` (seconds)

| Rep | P | chain_length input | outer | t_mcmc |
| --- | ---: | ---: | ---: | ---: |
| 1 | 100,000 | 223 | 1 | 76.515 |
| 1 | 100,000 | 446 | 2 | 84.821 |
| 1 | 200,000 | 223 | 1 | 147.150 |
| 1 | 200,000 | 446 | 2 | 171.600 |
| 2 | 100,000 | 223 | 1 | 70.930 |
| 2 | 100,000 | 446 | 2 | 84.973 |
| 2 | 200,000 | 223 | 1 | 142.255 |
| 2 | 200,000 | 446 | 2 | 169.454 |

### Extrapolated Target Time

- `SUMMARY_EST target_seconds = 3449.100`
- `target_hours = 0.958` (about 57.5 minutes)
- Replicate-based range: about `3274` to `3624` seconds (`0.91` to `1.01` hours)

## Standard BayesC Benchmark (`fast_blocks=false`)

### Design

- Slurm job: `30951596`
- Warmup run before timed cases
- Timed grid:
  - `P ∈ {100,000, 200,000}`
  - `chain_length ∈ {2, 4}`
  - `REPS = 2`
- Extrapolation:
  1. For each `P`, fit `t_mcmc = a + b * chain_length`
  2. Predict `t_mcmc` at `chain_length = 2000`
  3. Fit those predicted times vs `P` and extrapolate to `P = 2,000,000`

### Measured `t_mcmc` (seconds)

| Rep | P | chain_length | t_mcmc |
| --- | ---: | ---: | ---: |
| 1 | 100,000 | 2 | 27.457 |
| 1 | 100,000 | 4 | 40.214 |
| 1 | 200,000 | 2 | 55.186 |
| 1 | 200,000 | 4 | 107.764 |
| 2 | 100,000 | 2 | 35.555 |
| 2 | 100,000 | 4 | 53.252 |
| 2 | 200,000 | 2 | 56.058 |
| 2 | 200,000 | 4 | 109.762 |

### Extrapolated Target Time

- `SUMMARY_EST target_seconds = 735340.060`
- `target_hours = 204.261` (about 8.5 days)
- Replicate-based range: about `701539` to `769141` seconds (`194.9` to `213.7` hours)

## Summary Comparison

| Path | Estimated MCMC time at target (`N=50k`, `P=2M`, `L=2000`) |
| --- | --- |
| `fast_blocks=true` | `3449 sec` (`0.958 h`) |
| `fast_blocks=false` | `735340 sec` (`204.261 h`) |

Approximate speedup from block mode in this benchmark setup: about `213x`.

## Notes and Interpretation

- These are empirical estimates from short-chain calibration runs, not exact guarantees.
- For `fast_blocks=true`, chain interpretation is two-level (outer iterations + inner block repeats).
- The two Slurm jobs landed on different CPU node types, so absolute times include some hardware variance.
- Even with that caveat, the gap between block and non-block is very large for this target-scale scenario.
