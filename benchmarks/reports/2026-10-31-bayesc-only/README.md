# BayesC-only verification

Base: reworkhow/JWAS.jl master 038c6c8a497ec6e17e957657a5c13ae698377132,
JWAS 2.4.0. Runtime: Julia 1.11.7, Apple Silicon, one Julia and one BLAS thread.

## Scope and defaults

BayesC ordinary and fast-block updates support shared storage for one dense
genotype term in single- or multi-trait models. Both multi-trait samplers I/II
and sequential/independent blocks are covered. BayesR sampler files are unchanged.
Other methods keep their upstream input-copy/alignment/block-storage defaults.

BayesC defaults to matrix adoption and shared contiguous rows when the model
supports it. Input preprocessing can mutate matching-precision caller matrices.
Noncontiguous/repeated row selections need shared_genotypes=false; integer or
mismatched-precision matrix inputs need copy_genotypes=true for conversion.
File/DataFrame loading and unsupported model modes retain upstream behavior.

BayesC uses one in-place block residual path with multipliers matching the
genotype element type. There is no residual-mode option. Sequential MT updates
copy only the current alpha block and need no product workspace. The residual
buffer retains its element type. Other methods retain upstream arithmetic.
Rounding can alter seeded trajectories; the statistical model and schedule
are unchanged.

## Regression verification

All 1767 full-suite checks passed, covering the unchanged 988 legacy checks plus shared-memory,
ST/MT residual reconstruction, mixed-precision multipliers and legacy BayesB
arithmetic. Documentation built successfully; no deployment. Optional simulated-omics integration tests
and other runtime/OS combinations were not run.

## Direct prediction comparison

The same fixture and exporter were run against an untouched upstream archive and
the candidate, using the respective defaults. The 72 paired cases cover:

- single-trait BayesC and two-trait BayesC samplers I and II;
- both Float32 and Float64;
- ordinary, sequential and independent blocks;
- seeds 73, 123 and 456;
- annotated sequential cases for all trait/sampler/precision/seed combinations;
- one partially missing training trait and wholly held-out prediction animals.

Fixtures use 40 genotyped animals, 28 training animals, 19 markers, 500 outer
iterations, 150 burn-in, thinning 10, and 35 saved draws. The same priors,
fixed FARM and random LITTER effects, input precision and block schedule are used.

All 72 cases passed absolute/relative tolerance 1e-4 for EBV, PEV, marker
posterior mean/SD and model frequency. Maximum absolute differences were
1.108886e-6 (EBV), 8.178929e-7 (PEV), 5.1e-7 (marker mean), and 3.3e-7
(marker SD). Model frequencies matched exactly. All 18 ordinary cases matched
exactly; 20 of 54 block cases matched exactly. See summary.json and summary.csv.
Expected animal/marker counts and unique matching key sets were asserted before
joining by Trait/ID and Trait/Marker_ID. Marker quoting was normalized; sample
joined rows were inspected. Aligned tables are generated locally and ignored
in Git; summary.csv contains compact case-level evidence.

These synthetic finite-chain comparisons do not establish convergence,
prediction accuracy on production data, or universal bitwise equivalence.

## Storage comparison

Eight copied-versus-shared cases cover one/two traits, ordinary/sequential
blocks and two seeds. With 128 genotyped animals, 96 training animals, 257
markers, and output requested for all animals, retained genotype parent payload
fell from 230272 to 131584 bytes (42.9%) and parent count from two to one.
Outputs passed keyed comparisons. Both candidate paths use BayesC block views;
the comparison isolates row sharing. It excludes other model state and peak RSS.

## Reproduction

Run benchmarks/bayesc_prediction_compatibility.jl with output-directory and
upstream/candidate arguments, selecting the desired checkout with --project.
Compare exports with benchmarks/compare_bayesc_predictions.py followed by the
upstream directory, candidate directory and report directory.
Run benchmarks/bayesc_storage_comparison.jl for storage measurements.
The ordinary unit suite includes shared-memory and residual reconstruction checks.
Local execution used --compiled-modules=existing, --startup-file=no and
GKSwstype=100 to avoid cache writes and graphics initialization.
