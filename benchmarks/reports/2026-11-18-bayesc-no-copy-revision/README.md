# BayesC no-copy revision: production-path memory measurements

This report measures practical synthetic cases through the actual JWAS
`get_genotypes` → QC → model construction → alignment → single-trait block
BayesC MCMC path.
It does not establish feasibility or convergence for 50,000 individuals and
millions of SNPs.

## Measured results: 2026-11-18

All **36 production cases completed** and all **24 keyed comparisons passed**.
The exported EBV, PEV, marker mean, SD and inclusion frequencies matched exactly
between adopted/copying and adopted/file paths in every case (maximum exported
absolute difference zero). All 12 adopted cases retained their caller's matrix,
including every QC-on case. These are finite-chain regression results.

The run used Julia 1.11.7, Python 3.13.4, macOS 26.6.2 on Apple Silicon.
Production source content is recorded in local commit
`2de114ef501100cf5b4ab77ee6873a06e03fbe06`; the new benchmark scripts were
uncommitted during measurement. The recorded source-and-runner SHA-256 is
`05bb8345b24a74ffd79c08f1fb08f83ab2d872c8df70b7ffea462372e45dcb8d`.
Every case checked the exact worktree JWAS path and this content fingerprint.
The grid completed in 563.3 seconds with two workers; individual cold processes
took 28.2–36.7 seconds. Neither timing is a target-scale performance estimate.

The following are **absolute process high-water RSS values through MCMC**, in
MiB, showing the range across both seeds. They include startup and compilation.

| Individuals × SNPs | QC | Adopted peak MiB | Copied peak MiB | CSV peak MiB |
|---|---|---:|---:|---:|
| 128 × 257 | off | 915.6–943.4 | 860.0–914.9 | 842.2–872.0 |
| 128 × 257 | on | 873.2–909.6 | 889.3–926.8 | 819.1–823.4 |
| 512 × 1,025 | off | 940.5–1020.9 | 875.4–956.5 | 1011.6–1014.4 |
| 512 × 1,025 | on | 992.3–1034.1 | 901.2–982.0 | 890.4–980.6 |
| 2,048 × 4,097 | off | 993.7–1015.9 | 1038.9–1065.8 | 1014.5–1020.3 |
| 2,048 × 4,097 | on | 944.9–1007.3 | 956.2–1130.2 | 1061.8–1160.8 |

The package-load high-water baseline alone was 405.0–406.7 MiB. Across the grid,
post-import/QC high-water RSS ranged from 551.3–643.0 MiB for adoption,
542.5–665.3 MiB for explicit copying, and 556.5–748.9 MiB for CSV input.
The large startup/JIT contribution and run-to-run variation prevent a general
claim that adoption reduces total-process peak RSS by a fixed percentage.

Retained **genotype buffers**, unlike RSS, were deterministic across QC and seeds:

| Individuals × SNPs | Model parents: adopted/file | Model bytes: adopted/file | Model parents: copied | Model bytes: copied | Unique model + caller bytes: copied |
|---|---:|---:|---:|---:|---:|
| 128 × 257 | 1 | 131,584 | 2 | 230,272 | 361,856 |
| 512 × 1,025 | 1 | 2,099,200 | 2 | 3,673,600 | 5,772,800 |
| 2,048 × 4,097 | 1 | 33,562,624 | 2 | 58,734,592 | 92,297,216 |

The adopted model and caller share the same matrix, so their union equals the
adopted model bytes above. CSV input has no retained external caller matrix.
At the largest measured size, model genotype bytes fell 42.9% versus copying;
when the caller deliberately retains its original input, unique model-plus-input
bytes fell 63.6%. These percentages concern these buffer sets, not RSS. The
largest block Gram cache was 1,048,580 bytes in every path and is counted
separately. Full phase readings and timings are in `metrics.csv`.

## Reproduction and controls

From the repository root:

```sh
python3 benchmarks/run_bayesc_no_copy_rss.py --self-test
python3 benchmarks/run_bayesc_no_copy_rss.py
```

The default grid has 36 fresh Julia processes: 128 × 257, 512 × 1,025 and
2,048 × 4,097 genotype dimensions; matrix adoption, explicit matrix copying, and
the unchanged CSV file loader; QC off/on; and MCMC seeds 73 and 123. Each process
uses one Julia thread and one BLAS thread, `--startup-file=no`,
`--compiled-modules=existing`, 200 outer sweeps, 50 burn-in sweeps, thinning 10
(15 saved draws), and sequential 64-SNP blocks. The public input chain length is
200 × 64, following JWAS's existing block schedule. Two independent subprocesses
run concurrently; timings are descriptive and may reflect host contention.
`output_marker_effect_samples=false` suppresses per-draw marker-effect files;
posterior marker summaries and animal predictions are still produced. These
measurements therefore cover that production output setting, not all default
output files. `output_heritability=false`; the normal memory guard remains on.

The genotype fixture is reproducible and shared across seeds. Every SNP is
polymorphic and passes MAF 0.01, with no missing genotypes. QC-on therefore
exercises real imputation/MAF/variance code and the all-pass shortcut; it does
not benchmark actual SNP removal. Three quarters of the animals have
training-only phenotype records; predictions cover all animals. The MCMC
model and sampler arithmetic are identical across input cases. Two seeds
check reproducibility across storage paths, not convergence or accuracy.

Matrix cases deliberately retain the caller's original `X` through the end of
MCMC, verified by using that reference afterward. The adopted path shares row
and marker-block views; the copied comparator uses `copy_genotypes=true` and
`shared_genotypes=false`. Thus it compares both input ownership and row sharing.
The file case writes the same fixture before import, releases its fixture
matrix, calls GC, and imports with the existing CSV loader. Fixture generation
and file writing still contribute to the process-lifetime high-water marks.

The runner prints and checks the exact `pathof(JWAS)`. Metadata includes the
Git HEAD and a SHA-256 fingerprint of source, Project.toml and runner scripts,
so an uncommitted implementation is identifiable. Resumption requires matching
fingerprints and chain controls. Julia version and source path are also stored
in every metrics row.

## Interpreting the measurements

`metrics.csv` records absolute process-lifetime maximum RSS after package load,
before import, after import/QC, before MCMC, and after MCMC/alignment. Julia's
`Sys.maxrss()` is documented in the installed `base/sysinfo.jl` as **bytes**.
These are monotonic high-water marks, not current live RSS. They include Julia,
loaded packages, compilation, fixture preparation and runtime allocations.
Subtracting the package baseline does not yield exact genotype memory or an
isolated allocation count. Compilation can dominate these practical sizes;
smaller retained buffers do not guarantee smaller measured process peaks.

The buffer columns answer separate questions:

- `model_backing_bytes/count`: unique dense genotype parents retained by the
  JWAS genotype term, excluding other model state and external caller storage.
- `caller_matrix_bytes`: original matrix deliberately retained by the caller,
  or zero for the file case. This can alias model storage.
- `model_and_caller_unique_bytes`: the union of model genotype parents and
  caller input, counting aliases once. Block Gram caches are reported separately.
- `import_aliases_caller`: whether import/QC retains the exact matrix supplied
  by the caller. This is asserted for adoption with QC off and with all SNPs
  passing QC.

## Numerical comparison safeguards

Every animal and marker export first receives expected-count, unique-key and
expected-key-set checks. Marker keys normalize quotes and numeric SNP-index
forms through an explicit index-to-marker-label mapping. Predictions are joined
by `(Trait, ID)` and marker summaries by `(Trait, Marker_ID)`; validated joins
must preserve the expected row count. No row-order matching is used.

For each dimension/QC/seed combination, the adopted run is compared with copied
and file input. EBV, PEV, marker mean, SD and inclusion frequency use absolute
and relative tolerances of 1e-4. `comparisons.csv` contains maximum absolute
differences and pass flags. `joined_samples.csv` contains the first two keyed
joined rows for each comparison; complete aligned values and raw exports remain
under locally retained, gitignored `aligned/` and `raw/`. Logs are gitignored.
The comparator self-check confirms that reordered, quoted and numeric keys
align and that dropped, duplicated or incorrect keys are rejected.

## Scale interpretation and remaining allocation paths

At 50,000 individuals × 1,000,000 SNPs, a single Float32 genotype buffer is
200 GB decimal (186.3 GiB). Five million SNPs require 1 TB decimal for that
buffer alone. These are exact element-count calculations, not measured RSS or
runtime predictions. Float64 promotion needs a 400 GB destination per million
SNPs while the Float32 source remains alive, before other memory.

The revision preserves adoption when all SNPs pass QC. Actual SNP filtering
still allocates `genotypes[:, select]`; keeping nearly every SNP can transiently
require nearly two genotype buffers. The unchanged CSV loader materializes a
DataFrame, a dense Matrix and a mapped Matrix. Nonunit residual weights also
retain genotype-sized weighted columns. The normal MCMC memory guard executes
after import/QC/alignment, so it cannot protect those earlier peaks. The new
precision-conversion preflight checks known source/destination buffers before
promotion; it is not a total-process memory guarantee.

The measured grid is a regression and allocation-path check. A target-scale
feasibility claim requires a representative import and sufficiently sustained
production sampling run on the intended hardware, with peak RSS recorded.
